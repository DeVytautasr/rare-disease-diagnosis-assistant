#!/usr/bin/env python3
"""Phase 10: the attainable-ceiling echo on breakpoint_evidence_summary.

Phase 9 measured that no model at any tier could state the ceiling argument for
a heterozygous balanced rearrangement, because the score at which "strong"
begins appeared in no tool return and no tool description. This suite guards the
four properties that fix has to keep:

  1. the fields are present and derivable on a real call
  2. the band boundary is DERIVED from bam_tools' source, not a literal — if the
     scoring ladder is revised the reported boundary moves with it
  3. the change is PURELY ADDITIVE — no existing key removed or altered
  4. the echo fails soft and fails loud in the right places: it never breaks the
     tool it wraps, and an unparseable source is reported, never guessed

Run: python stage1_igv_assistant/tests/test_ceiling_echo.py
"""
import ast, asyncio, json, os, sys, textwrap

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(os.path.dirname(_HERE)))
from stage1_igv_assistant import score_tiers as st
from stage1_igv_assistant import server as ev

BAM = os.path.expanduser("~/public_data/sim/bams/IMP01.bam")
NEW_FIELDS = ["attainable_ceiling_derivable", "score_bands", "strong_band",
              "max_all_layers", "max_with_flat_depth", "attainable_here",
              "strong_band_reachable_here", "attainable_basis", "attainable_note"]
fails = []


def check(name, cond, detail=""):
    print(f"  {'PASS' if cond else 'FAIL'}  {name}{(' — ' + detail) if detail else ''}")
    if not cond:
        fails.append(name)


def summary(**kw):
    async def go():
        r = await ev.mcp.call_tool("breakpoint_evidence_summary",
                                   {"bam_path": BAM, "chromosome": "chr20",
                                    "position": 200000, **kw})
        return r.structured_content
    return asyncio.run(go())


print("1. fields present on a real call")
if not os.path.exists(BAM):
    print(f"  SKIP — {os.path.basename(BAM)} not present")
    d = None
else:
    d = summary()
    for f in NEW_FIELDS:
        check(f"returns {f}", f in d)
    check("derivable is True", d.get("attainable_ceiling_derivable") is True)
    check("strong_band is a number", isinstance(d.get("strong_band"), (int, float)))
    check("attainable_here <= max_all_layers",
          d["attainable_here"] <= d["max_all_layers"],
          f"{d['attainable_here']} <= {d['max_all_layers']}")
    check("strong_band_reachable_here agrees with the arithmetic",
          d["strong_band_reachable_here"] == (d["attainable_here"] >= d["strong_band"]))
    note = d.get("attainable_note", "")
    check("note states UNREACHABLE when the arithmetic says so",
          ("UNREACHABLE" in note) == (not d["strong_band_reachable_here"]), note[:60])
    check("note names the heterozygosity cap",
          "intact homolog" in note or "intact homolog" in d["attainable_basis"])
    check("basis names the top band's threshold, not a literal band name",
          str(int(max(t["threshold"] for t in
                      st.derive_tiers()["discordant_pairs"]["tiers"]))) in d["attainable_basis"]
          or "0.5" in d["attainable_basis"])

print("\n2. the band boundary is derived, not a literal")
bands = st.derive_bands()
check("derive_bands finds 'strong'", st.strong_band(bands) is not None,
      f"strong starts at {st.strong_band(bands)}")
if d:
    check("returned strong_band == derived strong_band",
          d["strong_band"] == st.strong_band(bands))
# a synthetic ladder with DIFFERENT numbers must yield those different numbers,
# which a hardcoded value could not do
fake = ast.parse(textwrap.dedent('''
    def summarize_breakpoint_evidence():
        if quality:
            evidence_strength = "QUALITY-LIMITED"
        elif evidence_score >= 88:
            evidence_strength = "strong"
        elif evidence_score >= 33:
            evidence_strength = "moderate"
        else:
            evidence_strength = "none"
'''))
fn = next(n for n in ast.walk(fake) if isinstance(n, ast.FunctionDef))
got = st._scan_bands(fn)
check("parser tracks a changed ladder (88/33, not 70/40)",
      got.get("strong", {}).get("threshold") == 88.0
      and got.get("moderate", {}).get("threshold") == 33.0,
      f"strong={got.get('strong', {}).get('threshold')} moderate={got.get('moderate', {}).get('threshold')}")
check("QUALITY-LIMITED is not treated as a band", "QUALITY-LIMITED" not in got)

print("\n3. purely additive")
if d:
    plain = ev.summarize_breakpoint_evidence(BAM, "chr20", 200000)
    wrapped = ev._with_ceiling(dict(plain))
    removed = sorted(set(plain) - set(wrapped))
    altered = sorted(k for k in set(plain) & set(wrapped)
                     if json.dumps(plain[k], sort_keys=True, default=str)
                     != json.dumps(wrapped[k], sort_keys=True, default=str))
    check("no key removed", not removed, str(removed))
    check("no key altered", not altered, str(altered))
    check("only the ceiling keys added", sorted(set(wrapped) - set(plain)) == sorted(NEW_FIELDS))

print("\n4. fails soft where it must, loud where it must")
check("error dicts pass through untouched",
      ev._with_ceiling({"error": "x", "error_type": "y"}) == {"error": "x", "error_type": "y"})
check("non-dicts pass through untouched", ev._with_ceiling(None) is None)
broken = ev._with_ceiling({"discordant_pairs": "not-a-dict", "soft_clips": None})
check("a malformed result does not raise", isinstance(broken, dict))
check("a malformed result still reports derivability",
      "attainable_ceiling_derivable" in broken)
try:
    empty = ast.parse("def summarize_breakpoint_evidence():\n    pass\n")
    fn2 = next(n for n in ast.walk(empty) if isinstance(n, ast.FunctionDef))
    check("a ladder-less source yields no bands", st._scan_bands(fn2) == {})
except Exception as e:
    check("a ladder-less source yields no bands", False, str(e))

print("\n" + "=" * 60)
print("ALL CEILING ECHO TESTS PASSED" if not fails else f"FAILED: {fails}")
sys.exit(1 if fails else 0)
