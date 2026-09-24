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

Sections 1 and 3 need a real BAM with a breakpoint (default: IMP01.bam at
chr20:200000). Without it they cannot run -- and a suite that tested nothing
must not report success. Every section counts the checks it executed; a section
that was skipped, a check that was skipped, or a section that ran zero checks
makes the run INCOMPLETE, which is reported loudly and has its own exit code.
This suite once printed "ALL CEILING ECHO TESTS PASSED" with section 1 skipped
and section 3 having run no assertions at all.

Exit codes:  0 = every section ran and every check passed
             1 = a check failed
             2 = INCOMPLETE: something was skipped or ran nothing (and nothing failed)

Data override (optional):  CEILING_ECHO_BAM=/path.bam  CEILING_ECHO_LOCUS=chr1:1050000

Run: python stage1_igv_assistant/tests/test_ceiling_echo.py
"""
import ast, asyncio, json, os, sys, textwrap

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(os.path.dirname(_HERE)))
from stage1_igv_assistant import score_tiers as st
from stage1_igv_assistant import server as ev

BAM = os.path.expanduser(os.environ.get("CEILING_ECHO_BAM", "~/public_data/sim/bams/IMP01.bam"))
CHROM, _, _POS = os.environ.get("CEILING_ECHO_LOCUS", "chr20:200000").rpartition(":")
POS = int(_POS.replace(",", ""))
NEW_FIELDS = ["attainable_ceiling_derivable", "score_bands", "strong_band",
              "max_all_layers", "max_with_flat_depth", "attainable_here",
              "strong_band_reachable_here", "attainable_basis", "attainable_note"]
fails = []
SECTIONS = []        # {"title", "checks", "skipped": reason or None, "skipped_checks": [...]}


def section(title):
    print(("\n" if SECTIONS else "") + title)
    SECTIONS.append({"title": title, "checks": 0, "skipped": None, "skipped_checks": []})


def skip(reason):
    """The whole current section cannot run."""
    SECTIONS[-1]["skipped"] = reason
    print(f"  SKIP  {reason}")


def skip_check(name, reason):
    """One check inside a section that otherwise runs."""
    SECTIONS[-1]["skipped_checks"].append(f"{name} ({reason})")
    print(f"  SKIP  {name} — {reason}")


def check(name, cond, detail=""):
    SECTIONS[-1]["checks"] += 1
    print(f"  {'PASS' if cond else 'FAIL'}  {name}{(' — ' + detail) if detail else ''}")
    if not cond:
        fails.append(name)


def summary(**kw):
    async def go():
        r = await ev.mcp.call_tool("breakpoint_evidence_summary",
                                   {"bam_path": BAM, "chromosome": CHROM,
                                    "position": POS, **kw})
        return r.structured_content
    return asyncio.run(go())


section("1. fields present on a real call")
if not os.path.exists(BAM):
    skip(f"{os.path.basename(BAM)} not present — no real call could be made")
    d = None
else:
    print(f"  using {os.path.basename(BAM)} at {CHROM}:{POS}")
    d = summary()
    for f in NEW_FIELDS:
        check(f"returns {f}", f in d)
    check("derivable is True", d.get("attainable_ceiling_derivable") is True)
    check("strong_band is a number", isinstance(d.get("strong_band"), (int, float)))
    check("attainable_here <= max_all_layers",
          d.get("attainable_here") is not None and d.get("max_all_layers") is not None
          and d["attainable_here"] <= d["max_all_layers"],
          f"{d.get('attainable_here')} <= {d.get('max_all_layers')}")
    check("strong_band_reachable_here agrees with the arithmetic",
          None not in (d.get("attainable_here"), d.get("strong_band"))
          and d.get("strong_band_reachable_here") == (d["attainable_here"] >= d["strong_band"]))
    note = d.get("attainable_note") or ""
    check("note states UNREACHABLE when the arithmetic says so",
          ("UNREACHABLE" in note) == (not d.get("strong_band_reachable_here")), note[:60])
    check("note names the heterozygosity cap",
          "intact homolog" in note or "intact homolog" in (d.get("attainable_basis") or ""))
    check("basis names the top band's threshold, not a literal band name",
          str(int(max(t["threshold"] for t in
                      st.derive_tiers()["discordant_pairs"]["tiers"]))) in (d.get("attainable_basis") or "")
          or "0.5" in (d.get("attainable_basis") or ""))

section("2. the band boundary is derived, not a literal")
bands = st.derive_bands()
check("derive_bands finds 'strong'", st.strong_band(bands) is not None,
      f"strong starts at {st.strong_band(bands)}")
if d:
    check("returned strong_band == derived strong_band",
          d.get("strong_band") == st.strong_band(bands))
else:
    skip_check("returned strong_band == derived strong_band", "needs section 1's real call")
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

section("3. purely additive")
if not d:
    skip("needs section 1's real call, which did not run")
else:
    plain = ev.summarize_breakpoint_evidence(BAM, CHROM, POS)
    wrapped = ev._with_ceiling(dict(plain))
    removed = sorted(set(plain) - set(wrapped))
    altered = sorted(k for k in set(plain) & set(wrapped)
                     if json.dumps(plain[k], sort_keys=True, default=str)
                     != json.dumps(wrapped[k], sort_keys=True, default=str))
    check("no key removed", not removed, str(removed))
    check("no key altered", not altered, str(altered))
    check("only the ceiling keys added", sorted(set(wrapped) - set(plain)) == sorted(NEW_FIELDS))

section("4. fails soft where it must, loud where it must")
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

# ── outcome: failed, incomplete, or complete — never "passed" by default ──
not_run = [s for s in SECTIONS if s["skipped"] or s["checks"] == 0]
partial = [s for s in SECTIONS if s["skipped_checks"] and s not in not_run]
ran = sum(s["checks"] for s in SECTIONS)
print("\n" + "=" * 60)
print(f"{ran} checks ran in {len(SECTIONS) - len(not_run)} of {len(SECTIONS)} sections")
if fails:
    print(f"FAILED: {fails}")
if not_run or partial:
    print("!" * 60)
    print("INCOMPLETE — this run did NOT test everything this suite covers:")
    for s in not_run:
        print(f"  NOT RUN  {s['title']}: {s['skipped'] or 'ran zero checks'}")
    for s in partial:
        for c in s["skipped_checks"]:
            print(f"  SKIPPED  {s['title']}: {c}")
    print("An incomplete run is not a pass. Provide the data (or set CEILING_ECHO_BAM /")
    print("CEILING_ECHO_LOCUS to a BAM with a breakpoint) and run it again.")
    print("!" * 60)
elif not fails:
    print("ALL CEILING ECHO TESTS PASSED")
sys.exit(1 if fails else (2 if (not_run or partial) else 0))
