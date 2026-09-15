#!/usr/bin/env python3
"""Phase 9 gate: nothing path-like or sample-like may reach the Anthropic API.

Every previous phase ran entirely on this machine, so a path in a tool return
was a cosmetic problem. It is not one any more. This test inspects the ACTUAL
serialised request body run_turn_api would send -- schemas, system prompt, user
message, and real recorded tool results -- and fails if any registered path,
any basename of one, the home directory, or the string 'patient_data' occurs
anywhere in it.

Both directions are controlled:
  POSITIVE  the payload built through _chat_exec is clean
  NEGATIVE  the same tool call with the scrub removed LEAKS, and the detector
            says so -- proving the positive result is not a detector that
            cannot fire.
Run: python stage1_igv_assistant/tests/test_api_leak.py
"""
import json, os, sys

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(_HERE))              # for stage1_igv_assistant.*
sys.path.insert(0, os.path.dirname(os.path.dirname(_HERE)))  # repo root
from stage1_igv_assistant import chat as chatmod
from stage1_igv_assistant import ui


def needles():
    """Every string whose appearance in an outbound payload is a leak."""
    n = set()
    for p in list(ui.DATASETS.values()) + list(ui.CANDIDATE_FILES.values()) + [ui.MASK_PATH]:
        n.add(p)
        n.add(os.path.basename(p))
        n.add(os.path.dirname(p))
    n.add(os.path.expanduser("~"))
    n.add("patient_data")
    n.add(os.path.expanduser("~/public_data"))
    return sorted((x for x in n if x and len(x) > 3), key=len, reverse=True)


def find_leaks(payload, needle_list):
    blob = json.dumps(payload, default=str)
    return sorted({x for x in needle_list if x in blob})


def main():
    ui.discover_public()
    assert ui.DATASETS, "no datasets registered; cannot test the real payload"
    N = needles()
    label = "NA12878.chr20_chr21" if "NA12878.chr20_chr21" in ui.DATASETS else sorted(ui.DATASETS)[0]
    print(f"registered: {len(ui.DATASETS)} datasets, {len(ui.CANDIDATE_FILES)} candidate files")
    print(f"leak needles: {len(N)} (longest {N[0]!r})")

    tools, where = ui._chat_tools(trim=False)
    api_tools = chatmod.to_anthropic_tools(tools)
    fails = 0

    # ---- 1. schemas -------------------------------------------------------
    hits = find_leaks(api_tools, N)
    print(f"\n[1] tool schemas ({len(api_tools)} tools) leaks={hits}")
    fails += bool(hits)
    enums = [p.get("enum") for t in api_tools for p in t["input_schema"].get("properties", {}).values()
             if isinstance(p, dict) and p.get("enum")]
    print(f"    label enums present: {len(enums)}; sample={enums[0][:3] if enums else None}")

    # ---- 2. a real tool call, through the real executor --------------------
    # applicable_layers echoes bam_path in its return -- it is the tool the
    # scrub exists for, so it is the one to test with.
    rec, err = ui._chat_exec("applicable_layers", {"dataset": label})
    assert err is None, f"executor refused: {err}"
    hits = find_leaks(rec["result"], N)
    print(f"\n[2] applicable_layers via _chat_exec  leaks={hits}")
    print(f"    result keys: {sorted(rec['result'])[:8]}")
    fails += bool(hits)

    # ---- 3. an ERROR return, which is where paths usually escape -----------
    bad, berr = ui._chat_exec("bam_stats_at_locus",
                              {"dataset": label, "chromosome": "chrNOPE",
                               "start": 1, "end": 500})
    body = bad["result"] if berr is None else {"error": berr}
    hits = find_leaks(body, N)
    print(f"\n[3] error return (bad contig)         leaks={hits}")
    print(f"    text: {json.dumps(body)[:180]}")
    fails += bool(hits)

    # ---- 4. THE WHOLE REQUEST BODY, as it would go on the wire -------------
    req = chatmod.run_turn_api(
        "claude-sonnet-5",
        f"Assess chr20:200000 in dataset {label} and report the evidence score.",
        tools, set(where), ui._chat_exec, dry_run=True)
    req = dict(req)
    req["messages"] = [{"role": "user", "content": req["messages"][0]["content"]},
                       {"role": "assistant", "content": "checking"},
                       {"role": "user", "content": [
                           {"type": "tool_result", "tool_use_id": "toolu_x",
                            "content": json.dumps(chatmod.shrink_for_model(
                                "applicable_layers", rec["result"]))},
                           {"type": "tool_result", "tool_use_id": "toolu_y",
                            "content": json.dumps(body)}]}]
    hits = find_leaks(req, N)
    size = len(json.dumps(req, default=str))
    print(f"\n[4] FULL REQUEST BODY ({size} bytes)   leaks={hits}")
    print(f"    keys: {sorted(req)}")
    print(f"    cache_control on system: "
          f"{req['system'][0].get('cache_control')}")
    fails += bool(hits)

    # ---- 5. NEGATIVE CONTROL: create the leak ------------------------------
    # Same tool, same coordinates, scrub bypassed. If this does NOT trip the
    # detector, tests 1-4 proved nothing.
    resolved, rerr = chatmod.resolve_args("applicable_layers", {"dataset": label},
                                          ui.DATASETS, ui.CANDIDATE_FILES, ui.MASK_PATH)
    assert rerr is None
    raw = ui.RECORDER.call(where["applicable_layers"], "applicable_layers", resolved)
    leaked_req = dict(req)
    leaked_req["messages"] = req["messages"][:2] + [{"role": "user", "content": [
        {"type": "tool_result", "tool_use_id": "toolu_z",
         "content": json.dumps(raw["result"])}]}]
    neg = find_leaks(leaked_req, N)
    print(f"\n[5] NEGATIVE CONTROL (scrub bypassed)  leaks={neg}")
    if not neg:
        print("    FAIL: the unscrubbed payload did not trip the detector.")
        fails += 1
    else:
        print(f"    detector fired on {len(neg)} needle(s) — it can fail. Good.")
        # and confirm the scrub is what removes it
        cleaned = find_leaks(ui.scrub(raw["result"]), N)
        print(f"    same result after ui.scrub(): leaks={cleaned}")
        fails += bool(cleaned)

    print("\nVERDICT:", "PASS — no leak on the API path" if fails == 0
          else f"FAIL ({fails} check(s) failed)")
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main())
