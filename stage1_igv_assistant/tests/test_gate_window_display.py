#!/usr/bin/env python3
"""
Regression test: the low-MAPQ figure the page shows must be measured over the
same window the quality gate decides on.

Named for the condition that exposed it. ui.assess() fetched the displayed
low_mapq_fraction over position +/-250 bp, while breakpoint_evidence_summary's
quality gate measures position +/-window_bp (500). At 21:19,281,000 in the
public HCC1143 slice the score was withheld with the tool reporting 43.6%,
while the page showed 35.6% -- a figure below the 40% gate, printed next to a
decision taken because the figure was above it. Six loci in that slice
straddled the gate this way, in both directions.

The fixture is synthetic and reproduces both directions:
  WITHHELD  low-MAPQ reads in the ring between 250 and 500 bp from the
            position: the +/-500 window is over the gate, +/-250 is not.
  SCORED    low-MAPQ reads in the core: +/-250 is over the gate, +/-500 is not.

Checks, per locus:
  CONDITION  the fixture really straddles: the +/-250 and +/-500 fractions,
             measured with the tool's own stats function, fall on opposite
             sides of the gate. If a fixture change ever stopped straddling,
             this test says so instead of passing vacuously.
  WINDOW     the window of the displayed figure equals the window the gate
             used, captured from INSIDE summarize_breakpoint_evidence -- so
             the test fails if either side's window moves without the other.
  AGREEMENT  the displayed figure is on the same side of the gate as the
             verdict, and the payload labels the window it was measured over.

Run: python3 stage1_igv_assistant/tests/test_gate_window_display.py
"""
import os
import sys
import tempfile

import pysam

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))

from stage1_igv_assistant import ui  # noqa: E402
from stage1_igv_assistant.tools import bam_tools  # noqa: E402
from stage1_igv_assistant.tools.bam_tools import LOW_MAPQ_QUALITY_GATE as GATE  # noqa: E402

FAILURES = []
HEADER = pysam.AlignmentHeader.from_dict({
    "HD": {"VN": "1.6", "SO": "coordinate"},
    "SQ": [{"SN": "chr1", "LN": 248956422}],
})
WINDOW_BP = 500      # what ui.assess passes to the summary tool by default

# locus name -> (position, MAPQ as a function of read start minus position)
LOCI = {
    "WITHHELD": (100_000, lambda d: 0 if (-600 < d <= -350 or 250 < d <= 500) else 60),
    "SCORED":   (200_000, lambda d: 0 if (-200 < d <= 100) else 60),
}


def check(label, condition, detail=""):
    if condition:
        print(f"  PASSED ✓  {label}")
    else:
        print(f"  FAILED ✗  {label}\n             {detail}")
        FAILURES.append(label)


def build(path):
    """Paired 100 bp reads every 5 bp across +/-3 kb of each locus."""
    n = 0
    with pysam.AlignmentFile(path, "wb", header=HEADER) as bam:
        for _, (p, mapq_of) in sorted(LOCI.items(), key=lambda kv: kv[1][0]):
            for start in range(p - 3000, p + 3000, 5):
                r = pysam.AlignedSegment(HEADER)
                r.query_name = f"r{n}"
                r.query_sequence = "A" * 100
                r.flag = 0x1 | 0x2
                r.reference_id = 0
                r.reference_start = start
                r.mapping_quality = mapq_of(start - p)
                r.cigar = [(0, 100)]
                r.next_reference_id = 0
                r.next_reference_start = start + 200
                r.template_length = 300
                r.query_qualities = pysam.qualitystring_to_array("I" * 100)
                bam.write(r)
                n += 1
    pysam.index(path)


def stats_fraction(bam, p, half):
    s = bam_tools.get_bam_stats_at_locus(bam, "chr1", max(0, p - half), p + half)
    return s["low_mapq_fraction"]


def assess_capturing_gate_window(label, p):
    """ui.assess, with every stats call made from inside the summary recorded."""
    captured = []
    real = bam_tools.get_bam_stats_at_locus

    def spy(bam_path, chromosome, start, end, *a, **k):
        captured.append((sys._getframe(1).f_code.co_name, chromosome, start, end))
        return real(bam_path, chromosome, start, end, *a, **k)

    bam_tools.get_bam_stats_at_locus = spy
    try:
        E = ui.assess(label, "chr1", p, window_bp=WINDOW_BP)
    finally:
        bam_tools.get_bam_stats_at_locus = real
    gate = [(s, e) for fn, _, s, e in captured if fn == "summarize_breakpoint_evidence"]
    return E, gate


def main():
    with tempfile.TemporaryDirectory() as tmp:
        bam = os.path.join(tmp, "gate_window.bam")
        build(bam)
        ui.DATASETS["GATEWIN"] = bam
        for name, (p, _) in LOCI.items():
            print(f"\n{name}: chr1:{p}")
            inner, outer = stats_fraction(bam, p, 250), stats_fraction(bam, p, WINDOW_BP)
            check(f"CONDITION fixture straddles the {GATE:.0%} gate "
                  f"(+/-250: {inner:.3f}, +/-{WINDOW_BP}: {outer:.3f})",
                  (inner > GATE) != (outer > GATE))

            E, gate = assess_capturing_gate_window("GATEWIN", p)
            if E.get("error"):
                check("assess returned a result", False, str(E.get("error")))
                continue
            S = E["summary"]
            shown = S.get("low_mapq_fraction")
            rec = ui.RECORDER.calls[E["stats"]["call"] - 1]
            shown_win = (rec["params"]["start"], rec["params"]["end"])
            check("WINDOW the gate was evaluated exactly once", len(gate) == 1, str(gate))
            check(f"WINDOW displayed figure's window {shown_win} == gate's window "
                  f"{gate[0] if gate else None}", gate and shown_win == gate[0])
            check("WINDOW payload labels the window it was measured over",
                  S.get("low_mapq_window") == {"start": shown_win[0], "end": shown_win[1],
                                               "half_width_bp": WINDOW_BP},
                  str(S.get("low_mapq_window")))
            check(f"AGREEMENT displayed figure {shown} is the gate-window figure {outer}",
                  shown == outer)
            withheld = S.get("evidence_strength") == "QUALITY-LIMITED"
            check(f"AGREEMENT displayed {shown} {'>' if shown > GATE else '<='} {GATE} "
                  f"and verdict is {S.get('evidence_strength')!r}",
                  (shown > GATE) == withheld)
            if name == "WITHHELD":
                check("the fixture drives the gate: score withheld",
                      withheld and S.get("evidence_score") is None, str(S.get("evidence_score")))
            else:
                check("the fixture scores normally: a number, not withheld",
                      not withheld and isinstance(S.get("evidence_score"), (int, float)),
                      f"{S.get('evidence_strength')!r} {S.get('evidence_score')!r}")

    print("\n" + "=" * 60)
    if FAILURES:
        print(f"{len(FAILURES)} FAILURE(S): {FAILURES}")
        sys.exit(1)
    print("ALL GATE-WINDOW DISPLAY TESTS PASSED")


if __name__ == "__main__":
    main()
