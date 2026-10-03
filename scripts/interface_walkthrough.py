#!/usr/bin/env python3
"""Click through the interface page in a real browser and record what it shows.

Phase 26: the page reviews every candidate at once (rearrangements sorted by the
reads that support them), shows no score, draws the rearrangement and the
abnormal reads at both ends, and names genes from the local table.

PUBLIC DATA ONLY. The private path is exercised with a PUBLIC file registered under
an explicit label (explicit registration is what makes a sample private), so the
consent gate, the refusal and the "sent with permission" line can all be tested
without any patient file being opened.

Start the interface first, from the repository root, and stop it by its PID after:

    python -m stage1_igv_assistant.ui \
        --dataset PRIVTEST=$HOME/public_data/sim/bams/IMP09.bam \
        --candidates PRIVTEST=$HOME/public_data/sim/delly/IMP09.bcf

then:

    python scripts/interface_walkthrough.py --out <dir> [--cloud]

Without --cloud nothing is sent to Anthropic: the cloud model is only selected, the
permission box is checked to appear, and one question that names private data is
posted WITHOUT permission, which the interface must refuse before any model runs.
With --cloud two questions are sent to claude-sonnet-5 (about 0.05-0.20 USD): one
about test data, and one about PRIVTEST with the permission box ticked.

Writes <out>/walkthrough.json (every fact below, with pass/fail) and numbered
screenshots. Exit status 1 if any check failed. Needs Playwright with Chromium
(pip install playwright; python -m playwright install chromium).
"""
import argparse
import asyncio
import json
import os
import sys
import time

from playwright.async_api import async_playwright


def parse():
    p = argparse.ArgumentParser()
    p.add_argument("--url", default="http://127.0.0.1:8765/")
    p.add_argument("--out", required=True)
    p.add_argument("--junction-text", default="14,100,001",
                   help="text identifying the candidate row to open (one end's position)")
    p.add_argument("--hand-dataset", default="IMP09")
    p.add_argument("--hand-position", default="chr20:33,700,000")
    p.add_argument("--private-label", default="PRIVTEST")
    p.add_argument("--private-position", default="chr20:33,700,000")
    p.add_argument("--local-model", default="qwen2.5:7b")
    p.add_argument("--cloud", action="store_true")
    p.add_argument("--no-igv", action="store_true")
    p.add_argument("--igv-timeout", type=int, default=900)
    p.add_argument("--width", type=int, default=1440)
    p.add_argument("--height", type=int, default=900)
    return p.parse_args()


async def main(a):
    os.makedirs(a.out, exist_ok=True)
    facts, checks, errors = {}, [], []
    shot_n = [0]

    def check(name, ok, detail=""):
        checks.append({"check": name, "ok": bool(ok), "detail": str(detail)[:300]})
        print(("  PASSED  " if ok else "  FAILED  ") + name + ("" if ok else f"  [{detail}]"), flush=True)

    async def shot(pg, name, full=False):
        shot_n[0] += 1
        await pg.screenshot(path=os.path.join(a.out, f"{shot_n[0]:02d}_{name}.png"), full_page=full)

    async def api(pg, method, path, body=None):
        return await pg.evaluate("""([m, p, b]) => fetch(p, b ? {method: m, headers: {'Content-Type': 'application/json'},
            body: JSON.stringify(b)} : {method: m}).then(r => r.json())""", [method, path, body])

    async with async_playwright() as p:
        b = await p.chromium.launch()
        pg = await b.new_page(viewport={"width": a.width, "height": a.height})
        pg.on("pageerror", lambda e: errors.append("pageerror: " + str(e)))
        pg.on("console", lambda m: errors.append("console: " + m.text) if m.type == "error" else None)

        # 1. the page opens on test data, and loads nothing private by itself
        print("1. candidates", flush=True)
        await pg.goto(a.url)
        await pg.wait_for_selector("#funnel .funnel-count", timeout=120000)
        await pg.wait_for_selector("tr.rearr, #cands .err", timeout=300000)
        await pg.wait_for_timeout(800)
        facts["default_sample"] = await pg.input_value("#sample")
        facts["default_sample_kind"] = await pg.inner_text("#sampleKind")
        check("the page opens on a test-data sample", facts["default_sample_kind"] == "Test data", facts)
        boot = await api(pg, "GET", "/api/bootstrap")
        facts["kinds"] = boot.get("kinds")
        check(f"{a.private_label} (registered explicitly) is reported as private",
              (boot.get("kinds") or {}).get("datasets", {}).get(a.private_label) == "private", boot.get("kinds"))
        calls = (await api(pg, "GET", "/api/calls")).get("calls", [])
        loaded = [c["params"].get("label") for c in calls if c["tool"] == "load_candidate_set"]
        facts["loaded_at_start"] = loaded
        check("nothing private was loaded at start", a.private_label not in loaded, loaded)
        facts["funnel_counts"] = [t.strip().split("\n")[0] for t in await pg.locator("#funnel .funnel-count").all_inner_texts()]
        facts["candidates_title"] = await pg.inner_text("#cands-title")
        facts["rearrangement_rows"] = await pg.locator("tr.rearr").count()
        facts["first_rows"] = [" | ".join(t.split()) for t in (await pg.locator("tr.rearr").all_inner_texts())[:5]]
        facts["funnel_explanation"] = await pg.inner_text("#funnel-why")
        check("every candidate was reviewed: rearrangements are listed, most read support first",
              facts["rearrangement_rows"] > 0, facts["candidates_title"])
        view = (await pg.inner_text("#view-cands")).lower()
        check("the candidate list shows no score", not any(w in view for w in ("score of", "combined score", "moderate 40")),
              [w for w in ("score of", "combined score", "moderate 40") if w in view])
        await shot(pg, "candidates", full=True)

        # 2. one junction's evidence
        print("2. evidence", flush=True)
        row = pg.locator("tr.rearr", has_text=a.junction_text).first
        check(f"a rearrangement containing {a.junction_text} is listed", await row.count() == 1)
        facts["row"] = " | ".join((await row.inner_text()).split())
        await row.locator("button", has_text="Review").click()
        await pg.wait_for_selector(".readbox svg, #reads-panel .err", timeout=300000)
        await pg.wait_for_timeout(1500)
        facts["evidence_title"] = await pg.inner_text(".ev-title")
        facts["summary_line"] = await pg.inner_text(".summary-line")
        facts["cautions"] = await pg.locator(".cautions li").all_inner_texts()
        facts["genes"] = " ".join((await pg.inner_text(".gtable") if await pg.locator(".gtable").count()
                                   else await pg.locator(".panel", has_text="Genes").first.inner_text()).split())[:600]
        facts["read_box_heads"] = await pg.locator(".readbox h4").all_inner_texts()
        facts["legend"] = await pg.locator(".legend").first.inner_text() if await pg.locator(".legend").count() else None
        check("the rearrangement is drawn: the chromosomes before and after", await pg.locator(".schematic svg").count() == 1)
        check("the abnormal reads at both ends are drawn", await pg.locator(".readbox svg").count() >= 2)
        ev = (await pg.inner_text("#view-evidence")).lower()
        check("the evidence view shows no score", not any(w in ev for w in ("score", "reachable here", "moderate 40", "strong 70")),
              [w for w in ("score", "reachable here", "moderate 40", "strong 70") if w in ev])
        await shot(pg, "evidence", full=True)

        # 3. IGV (Phase 27): one image of both ends, with only the reads the page draws
        if not a.no_igv:
            print("3. IGV (a minute)", flush=True)
            t0 = time.time()
            facts["igv_buttons"] = await pg.locator("[data-igv]").all_inner_texts()
            await pg.locator('[data-igv="0"]').click()
            await pg.wait_for_timeout(800)
            check("while IGV draws, its button is disabled", await pg.locator('[data-igv="0"]').is_disabled())
            await shot(pg, "igv_waiting")
            await pg.wait_for_selector("#igvb-0 .igv-grid, #igvb-0 .err", timeout=a.igv_timeout * 1000)
            await pg.wait_for_timeout(1500)
            facts["igv_seconds"] = round(time.time() - t0)
            facts["igv_images"] = await pg.locator("#igvb-0 .igv-grid img").count()
            facts["igv_caption"] = " ".join((await pg.inner_text("#igvb-0")).split())[:600]
            calls = (await api(pg, "GET", "/api/calls")).get("calls", [])
            view = next((c for c in reversed(calls) if c["tool"] == "igv_junction_view"), None)
            res = ((await api(pg, "GET", f"/api/call?id={view['id']}")) if view else {}).get("result") or {}
            facts["igv_result"] = {k: res.get(k) for k in ("windows", "reads_drawn", "reads_not_found_again",
                                                           "groups", "normal_reads_hidden", "genes_shown")}
            check("IGV made one image of both ends", facts["igv_images"] == 1 and len(res.get("windows") or []) == 2,
                  await pg.inner_text("#igvb-0"))
            check("... with only the reads the page draws, none of them missing",
                  res.get("reads_drawn", 0) > 0 and res.get("reads_not_found_again") == 0, facts["igv_result"])
            await shot(pg, "igv", full=True)

        # 4. a position typed in by hand
        print("4. hand entry", flush=True)
        await pg.select_option("#cp-ds", a.hand_dataset)
        await pg.fill("#cp-pos", a.hand_position)
        await pg.click("#checkpos button")
        await pg.wait_for_selector(".banner.hand, #ev .err", timeout=120000)
        await pg.wait_for_timeout(1000)
        facts["hand_title"] = await pg.inner_text(".ev-title")
        facts["hand_partners"] = " ".join((await pg.inner_text("#ev")).split())[:500]
        check("a position typed in by hand says where its reads point, without a score",
              "Where the reads here point" in await pg.inner_text("#ev") and "score" not in (await pg.inner_text("#ev")).lower())
        await shot(pg, "hand_entry", full=True)
        # the same position in the private-labelled copy, so it becomes a private position
        await pg.select_option("#cp-ds", a.private_label)
        await pg.fill("#cp-pos", a.private_position)
        await pg.click("#checkpos button")
        await pg.wait_for_function("lbl => document.querySelector('#ev-head').innerText.includes(lbl) && "
                                   "document.querySelector('.banner.hand, #ev .err')", arg=a.private_label, timeout=120000)
        await pg.wait_for_timeout(800)
        check("a private sample is marked private next to its name",
              "Private data" in await pg.inner_text("#ev-head"))

        # 5. the assistant on this computer
        print("5. assistant, local", flush=True)
        await pg.click("#tab-assistant")
        await pg.wait_for_timeout(600)
        models = await api(pg, "GET", "/api/chat_models")
        facts["local_models"], facts["cloud_models"] = models.get("models"), models.get("cloud_models")
        if a.local_model in (models.get("models") or []):
            await pg.select_option("#model", a.local_model)
            await pg.fill("#question", f"Which rearrangements in {facts['default_sample']} have the most read support, "
                                       f"and which genes do they break?")
            await pg.click("#askbtn")
            await pg.wait_for_selector(".qa .answer, .qa .err, .qa .verdict", timeout=900000)
            await pg.wait_for_timeout(800)
            facts["local_verdict"] = (await pg.locator(".qa .verdict").all_inner_texts())[:3]
            facts["local_tool_calls"] = await pg.locator(".qa .call").count()
            facts["local_tool_names"] = await pg.locator(".qa .call .name").all_inner_texts()
            facts["local_answer"] = (await pg.locator(".qa .answer").first.inner_text())[:800] if await pg.locator(".qa .answer").count() else None
            check("the local model answered through tool calls", facts["local_tool_calls"] >= 1, facts["local_verdict"])
            check("... and reviewed every candidate in one call", "review_candidates" in facts["local_tool_names"],
                  facts["local_tool_names"])
            await shot(pg, "assistant_local")
        else:
            check(f"local model {a.local_model} is available", False, models.get("why"))

        # 6. the cloud permission gate (nothing is sent here)
        print("6. permission gate", flush=True)
        if "claude-sonnet-5" in (models.get("cloud_models") or []):
            await pg.select_option("#model", "claude-sonnet-5")
            await pg.fill("#question", f"In dataset {a.private_label}, where do the reads at {a.private_position} point?")
            await pg.dispatch_event("#question", "input")
            await pg.wait_for_timeout(900)
            check("naming a private sample shows the permission box", await pg.locator("#consent").count() == 1)
            check("... and the Ask button stays off until it is ticked", await pg.locator("#askbtn").is_disabled())
            await shot(pg, "permission_named")
            pos_only = f"Which gene, if any, lies at {a.private_position.replace(',', '')}?"
            await pg.fill("#question", pos_only)
            await pg.dispatch_event("#question", "input")
            await pg.wait_for_timeout(900)
            check("a position read from private data also needs permission", await pg.locator("#consent").count() == 1,
                  await pg.inner_text("#privacy"))
            r = await api(pg, "POST", "/api/chat", {"model": "claude-sonnet-5", "message": pos_only})
            facts["refusal"] = (r.get("error") or "")[:200]
            check("posted without permission, it is refused before any model runs",
                  facts["refusal"].startswith("Not sent") and not r.get("events"), facts["refusal"])
            await pg.fill("#question", f"Which rearrangements in {facts['default_sample']} have the most read support, "
                                       f"and which genes do they break?")
            await pg.dispatch_event("#question", "input")
            await pg.wait_for_timeout(900)
            check("a test-data question needs no permission",
                  await pg.locator("#consent").count() == 0 and not await pg.locator("#askbtn").is_disabled())
        else:
            check("claude-sonnet-5 is offered (needs the API key)", False, models.get("cloud_models"))

        # 7. optional: two real cloud answers, test data only
        if a.cloud and "claude-sonnet-5" in (models.get("cloud_models") or []):
            print("7. assistant, cloud (sends two questions)", flush=True)
            n0 = await pg.locator(".qa .meta").count()
            await pg.click("#askbtn")
            await pg.wait_for_function(f"document.querySelectorAll('.qa .meta').length > {n0}", timeout=900000)
            await pg.wait_for_timeout(800)
            first = pg.locator(".qa").first
            facts["cloud_public_meta"] = await first.locator(".meta").inner_text()
            facts["cloud_public_lines"] = await first.locator(".verdict").all_inner_texts()
            check("the cloud model answered a test-data question", await first.locator(".answer").count() == 1,
                  facts["cloud_public_lines"])
            await shot(pg, "assistant_cloud")
            n = await pg.locator(".qa .meta").count()
            await pg.fill("#question", f"In dataset {a.private_label}, where do the reads at {a.private_position} point?")
            await pg.dispatch_event("#question", "input")
            await pg.wait_for_timeout(900)
            await pg.check("#consent")
            await pg.click("#askbtn")
            await pg.wait_for_function(f"document.querySelectorAll('.qa .meta').length > {n}", timeout=900000)
            await pg.wait_for_timeout(800)
            first = pg.locator(".qa").first
            facts["cloud_private_meta"] = await first.locator(".meta").inner_text()   # time, turns, tokens, cost
            facts["cloud_private_lines"] = await first.locator(".verdict").all_inner_texts()
            check("with permission, the answer says what was sent",
                  any("Sent to Anthropic with your permission" in t for t in facts["cloud_private_lines"]),
                  facts["cloud_private_lines"])
            # permission is per question: the same text is still in the box, so the
            # permission box is still shown, and it must be unticked again
            await pg.wait_for_selector("#consent", timeout=20000)
            facts["consent_ticked_after_sending"] = await pg.locator("#consent").is_checked()
            check("once the question is sent, the permission box is unticked again",
                  facts["consent_ticked_after_sending"] is False and await pg.locator("#askbtn").is_disabled())
            await shot(pg, "assistant_cloud_permission")

        # 8. a private sample's list stays hidden until asked for
        print("8. private sample list", flush=True)
        await pg.click("#tab-cands")
        n_calls = len((await api(pg, "GET", "/api/calls")).get("calls", []))
        await pg.select_option("#sample", a.private_label)
        await pg.wait_for_selector("#reveal", timeout=120000)
        await pg.wait_for_timeout(600)
        facts["private_hidden_title"] = await pg.inner_text("#cands-title")
        check("choosing a private sample shows its counts but no candidate row",
              await pg.locator("tr.cand, tr.rearr").count() == 0 and await pg.locator("#funnel .funnel-count").count() > 0)
        calls = (await api(pg, "GET", "/api/calls")).get("calls", [])[n_calls:]
        lists = [c for c in calls if c["tool"] == "list_candidates"]
        check("... and the private sample is not reviewed before it is shown",
              not [c for c in calls if c["tool"] == "review_candidates"], [c["tool"] for c in calls])
        check("... and asks the server for counts only (no position reaches the page or the call log)",
              bool(lists) and all((c.get("params") or {}).get("offset", 0) >= 2 ** 31 - 1 for c in lists),
              [(c.get("params") or {}).get("offset") for c in lists])
        await shot(pg, "private_list_hidden", full=True)
        # a comparison chosen while the list is still hidden, as the demonstration does:
        # the count changes, no row appears, and the reply to the page carries no position
        # (the server's two listing calls for the comparison are in the call log)
        await pg.select_option("#f-rec", facts["default_sample"])
        await pg.wait_for_function("document.querySelector('#funnel').innerText.includes('Not found in')", timeout=120000)
        await pg.wait_for_timeout(600)
        held = await pg.evaluate("JSON.stringify(S.compare)")
        facts["private_hidden_with_comparison"] = {
            "title": await pg.inner_text("#cands-title"),
            "button": await pg.inner_text("#reveal") if await pg.locator("#reveal").count() else None,
            "last_funnel_row": " ".join((await pg.locator("#funnel .step").last.inner_text()).split())}
        check("with a comparison chosen the list stays hidden, and the page is sent no position",
              await pg.locator("tr.cand, tr.rearr").count() == 0 and await pg.locator("#reveal").count() == 1
              and '"survivor_effect"' in held and '"pos1"' not in held, held[:300])
        await pg.select_option("#f-rec", "")
        await pg.wait_for_function("!document.querySelector('#funnel').innerText.includes('Not found in')", timeout=120000)
        await pg.wait_for_timeout(600)
        # with the page's defaults (translocations only) the private copy can have no
        # junction at all (IMP09 has none): the reveal then shows that nothing passes
        facts["private_matching"] = await pg.evaluate("S.funnel.total_matching")
        await pg.click("#reveal")
        await pg.wait_for_selector("tr.rearr, #cands .err" + (", #cands .empty" if facts["private_matching"] == 0 else ""),
                                   timeout=300000)
        await pg.wait_for_timeout(600)
        if facts["private_matching"]:
            check("Show the candidates reveals the list, reviewed", await pg.locator("tr.rearr").count() > 0)
        else:
            facts["private_revealed"] = await pg.inner_text("#cands")
            check("Show the candidates reveals the list (empty here: nothing passes the filters)",
                  "No junction passes these filters" in facts["private_revealed"], facts["private_revealed"])

        # 9. how an answer is drawn: the page's own function on a text written here (no model)
        print("9. answer rendering", flush=True)
        drawn = await pg.evaluate("""(t) => { const d = document.createElement('div'); d.className = 'answer';
            d.innerHTML = renderMarkdown(t, [{text: '39', supported: true}]); document.body.appendChild(d);
            const numbers = [...d.querySelectorAll('ol > li')].map(li => { const ol = li.parentElement;
                return li.hasAttribute('value') ? li.value : (ol.start || 1) + [...ol.children].indexOf(li); });
            const out = {numbers_shown: numbers, text: d.innerText,
                         elements_from_model_html: d.querySelectorAll('img,script,a,u,iframe').length};
            d.remove(); return out; }""",
            "Two steps:\n\n1. It's 39 pairs\n\n2. <img src=x onerror=alert(1)> <u>second</u>\n\n3. third")
        facts["answer_rendering"] = drawn
        check("a numbered list keeps the numbers the model wrote, also across blank lines",
              drawn["numbers_shown"] == [1, 2, 3], drawn)
        check("an apostrophe next to a checked 39 stays an apostrophe, and HTML the model wrote stays text",
              "It's 39 pairs" in drawn["text"] and "<img" in drawn["text"] and drawn["elements_from_model_html"] == 0, drawn)

        # 10. limits and the call log
        await pg.evaluate("window.scrollTo(0,0)")
        await pg.click("#btnLimits")
        await pg.wait_for_timeout(500)
        facts["limits_headings"] = await pg.locator("#modal .limit h3").all_inner_texts()
        await shot(pg, "limits")
        await pg.keyboard.press("Escape")
        await pg.click("#btnLog")
        await pg.wait_for_timeout(600)
        await shot(pg, "call_log")
        await pg.keyboard.press("Escape")
        await b.close()

    check("no script error on the page", not errors, errors[:3])
    out = {"url": a.url, "cloud_questions_sent": bool(a.cloud), "facts": facts, "checks": checks, "page_errors": errors}
    with open(os.path.join(a.out, "walkthrough.json"), "w") as f:
        json.dump(out, f, indent=1)
    failed = [c["check"] for c in checks if not c["ok"]]
    print(f"\n{len(checks) - len(failed)} of {len(checks)} checks passed" + (f"; FAILED: {failed}" if failed else ""))
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(asyncio.run(main(parse())))
