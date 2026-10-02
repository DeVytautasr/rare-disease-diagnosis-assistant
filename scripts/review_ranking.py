#!/usr/bin/env python3
"""Phase 26: where the implanted translocations land in the review's ranking.

PUBLIC DATA ONLY: the twelve synthetic samples of the controlled positive test
(results/synthetic_control_2026-09/implants_ground_truth.json) and the background
they were built on.

    python3 scripts/review_ranking.py OUT_JSON [--sim ~/public_data/sim]
           [--background-bam ~/public_data/NA12878.chr20_chr21.bam]
           [--background-calls ~/public_data/sim/delly/background.bcf]

For every sample with both a read file (sim/bams/IMPxx.bam) and a caller file
(sim/delly/IMPxx.bcf) it runs review_candidates exactly as the page does (its
defaults: translocations, PASS, read pairs >= 3, split reads >= 1, main
chromosomes, the configured exclude template), and records:

  detected        a rearrangement has a breakpoint within MATCH_BP of the
                  implant's chr20 breakpoint and one within MATCH_BP of its chr21
                  breakpoint
  rank            its place in the list (1 = most supporting reads), and the
                  number of rearrangements listed
  pattern         "both junctions" when both of the implant's junctions passed
                  the filters, otherwise what the review says
  support         read pairs, split reads and distinct reads it counted
  cautions        the cautions the page would show

The background (no implant) gives the list a geneticist would see in a sample
without the event: its length and its best-supported rearrangement.
Nothing is written into the repository except OUT_JSON, which holds public
coordinates only.
"""
import argparse
import json
import os
import sys
import time

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, REPO)

from stage1_igv_assistant import config as CFG                    # noqa: E402
from stage1_igv_assistant import review_server as R               # noqa: E402
from stage1_igv_assistant.tools import vcf_tools                  # noqa: E402

TRUTH = os.path.join(REPO, "stage1_igv_assistant", "results", "synthetic_control_2026-09",
                     "implants_ground_truth.json")
MATCH_BP = 1000      # author judgement: the review's own read window


def tilde(p):
    h = os.path.expanduser("~")
    return p.replace(h, "~", 1) if p and p.startswith(h) else p


def matches(ev, bps):
    def near(chrom, pos):
        return any(b["chromosome"] == chrom and b["start"] - MATCH_BP <= pos <= b["end"] + MATCH_BP
                   for b in ev["breakpoints"])
    return all(near(c, p) for c, p in bps.items())


def summarise(ev):
    return {"rank": ev["rank"], "pattern": ev["pattern"], "iscn": ev.get("iscn"),
            "breakpoints": [{k: b[k] for k in ("chromosome", "start", "end", "band")} for b in ev["breakpoints"]],
            "support": ev["support"], "cautions": ev.get("cautions", []),
            "junctions": [{"caller_id": j.get("caller_id"), "orientation": j.get("orientation"),
                           "pe": j.get("pe"), "sr": j.get("sr"), "support": j.get("support")}
                          for j in ev["junctions"]]}


def review(calls, bam):
    vcf_tools.reset_registry()
    t0 = time.time()
    r = R.review_candidates(calls, bam)
    return r, round(time.time() - t0, 2)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("out")
    ap.add_argument("--sim", default="~/public_data/sim")
    ap.add_argument("--background-bam", default="~/public_data/NA12878.chr20_chr21.bam")
    ap.add_argument("--background-calls", default="~/public_data/sim/delly/background.bcf")
    a = ap.parse_args()
    if os.path.exists(a.out):
        sys.exit("the record exists")
    mask = CFG.status("exclude_template", kind="file")
    genes = CFG.status("gene_table", kind="file")
    R.configure(mask_path=mask["path"] if mask["found"] else None,
                gene_table=genes["path"] if genes["found"] else None,
                mim2gene=CFG.status("omim_mim2gene")["path"], gene_disorders=CFG.status("gene_disorders")["path"])
    truth = json.load(open(TRUTH))
    sim = os.path.expanduser(a.sim)
    rec = {"what": "rank of each implanted translocation in the review's list (Phase 26)",
           "definitions": " ".join(__doc__.split()), "started": time.strftime("%Y-%m-%dT%H:%M:%S"),
           "match_bp": MATCH_BP, "exclude_template_found": mask["found"], "gene_table_found": genes["found"],
           "samples": [], "background": None}
    for imp in truth["implants"]:
        sid = imp["id"]
        bam = os.path.join(sim, "bams", f"{sid}.bam")
        calls = os.path.join(sim, "delly", f"{sid}.bcf")
        row = {"id": sid, "class": imp.get("class"), "breakpoints": imp["breakpoints"]}
        if not (os.path.isfile(bam) and os.path.isfile(calls)):
            row["not_run"] = f"missing {tilde(bam) if not os.path.isfile(bam) else tilde(calls)}"
            rec["samples"].append(row)
            print(sid, row["not_run"], flush=True)
            continue
        r, secs = review(calls, bam)
        if "error" in r:
            row["error"] = r["error"]
            rec["samples"].append(row)
            continue
        bps = {"chr20": imp["breakpoints"]["chr20"], "chr21": imp["breakpoints"]["chr21"]}
        hit = next((ev for ev in r["rearrangements"] if matches(ev, bps)), None)
        row.update({"seconds": secs, "junctions_passing_filters": r["junctions_passing_filters"],
                    "rearrangements_listed": r["rearrangements_found"], "detected": hit is not None,
                    "implant": summarise(hit) if hit else None,
                    "first": summarise(r["rearrangements"][0]) if r["rearrangements"] else None,
                    "patterns": {p: sum(1 for e in r["rearrangements"] if e["pattern"] == p)
                                 for p in sorted({e["pattern"] for e in r["rearrangements"]})}})
        rec["samples"].append(row)
        print(sid, "detected" if hit else "not detected", f"rank {hit['rank']} of {r['rearrangements_found']}"
              if hit else f"({r['rearrangements_found']} listed)", hit and hit["pattern"], flush=True)
    bg_bam, bg_calls = os.path.expanduser(a.background_bam), os.path.expanduser(a.background_calls)
    if os.path.isfile(bg_bam) and os.path.isfile(bg_calls):
        r, secs = review(bg_calls, bg_bam)
        rec["background"] = {"seconds": secs, "junctions_passing_filters": r.get("junctions_passing_filters"),
                             "rearrangements_listed": r.get("rearrangements_found"),
                             "first": summarise(r["rearrangements"][0]) if r.get("rearrangements") else None,
                             "supports": [e["support"]["fragments"] for e in r.get("rearrangements", [])][:20]}
    else:
        rec["background"] = {"not_run": "background read or caller file missing"}
    det = [s for s in rec["samples"] if s.get("detected")]
    rec["summary"] = {
        "samples_run": sum(1 for s in rec["samples"] if "rearrangements_listed" in s),
        "implants_detected": len(det),
        "ranked_first": sum(1 for s in det if s["implant"]["rank"] == 1),
        "ranks": {s["id"]: s["implant"]["rank"] for s in det},
        "both_junctions": sum(1 for s in det if s["implant"]["pattern"] == "both junctions"),
        "background_listed": (rec["background"] or {}).get("rearrangements_listed")}
    rec["finished"] = time.strftime("%Y-%m-%dT%H:%M:%S")
    json.dump(rec, open(a.out, "w"), indent=1)
    print(json.dumps(rec["summary"]))


if __name__ == "__main__":
    main()
