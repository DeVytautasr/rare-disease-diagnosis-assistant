#!/usr/bin/env python3
"""Build a demo bundle small enough to send.

The whole-chromosome BAMs are 1.6 GB each, which is the only thing making a
demo impractical -- the candidate BCFs are 127 KB. So the BAMs are sliced to a
few windows around known breakpoints. Slicing is done with pysam, not samtools:
this tool never invokes samtools and the bundle must not require it either.

WHAT A SLICE PRESERVES
  * the full BAM header, so contig validation still behaves exactly as before
  * every read overlapping the kept windows, so all four evidence layers, the
    reciprocal check and the scoring behave exactly as on the full file
  * the candidate BCFs, untouched and whole -- the filter chain and
    compare_candidate_sets read only those and are unaffected

WHAT A SLICE BREAKS
  * any locus OUTSIDE a kept window returns zero reads. It does not error: the
    layers report assessable=false, which is honest but is not the same picture
    the full BAM gives. Only the loci listed in DEMO.md are meaningful.
  * genome-wide depth statistics are meaningless on a slice.
"""
import os, sys, pysam

WINDOW = 10_000

BUNDLE = {
    "DEMO_CLEAN": {
        "bam": "~/public_data/sim/bams/IMP01.bam",
        "bcf": "~/public_data/sim/delly/IMP01.bcf",
        "regions": [("chr20", 200_000), ("chr21", 14_100_000)],
        "why": "clean, uniquely mappable breakpoint",
    },
    "DEMO_REPEAT": {
        "bam": "~/public_data/sim/bams/IMP10.bam",
        "bcf": "~/public_data/sim/delly/IMP10.bcf",
        "regions": [("chr20", 25_800_000), ("chr21", 7_600_000)],
        "why": "low-mappability breakpoint",
    },
}


def measured(bam, chrom, pos):
    """What DEMO.md says about a locus, measured on the slice just written.

    DEMO.md used to carry literals ("all four layers; score 47.5/100") from the
    first implant run; when the implants were rebuilt they no longer described the
    data the bundle held. Now the numbers come from breakpoint_evidence_summary on
    the sliced BAM itself, through the MCP dispatch and the applicable-layer
    sampling (20,000 reads) the interface uses, so the text cannot drift from the
    file beside it."""
    import asyncio
    from stage1_igv_assistant import server

    def call(tool, params):
        return asyncio.run(server.mcp.call_tool(tool, params)).structured_content
    layers = call("applicable_layers", {"bam_path": bam, "sample_reads": 20000}).get("applicable_layers")
    s = call("breakpoint_evidence_summary", {"bam_path": bam, "chromosome": chrom, "position": pos,
                                             "label": f"{chrom}:{pos}", "applicable_layers": layers,
                                             "window_bp": 500})
    if "error" in s:
        return f"the tool returned an error: {s['error']}"
    comps = [s.get(k) for k in ("discordant_pair_score", "soft_clip_score", "split_read_score", "depth_score")]
    fired = sum(1 for c in comps if c)
    if s["evidence_strength"] == "QUALITY-LIMITED":
        return (f"score WITHHELD, `QUALITY-LIMITED` — not a low score ({s['locus_stats']['low_mapq_fraction']:.0%} "
                f"of reads below MAPQ 20); {fired} of 4 layers show signal")
    out = f"{fired} of 4 layers show signal; score {s['evidence_score']}/100 `{s['evidence_strength']}`"
    if s.get("strong_band_reachable_here") is False:
        out += (f"; the top band is unreachable for a balanced event (at most {s['attainable_here']:g} "
                f"reachable here, {s['strong_band']:g} needed)")
    return out


def slice_bam(src, dst, regions, window=WINDOW):
    src = os.path.expanduser(src)
    n_in = n_out = 0
    with pysam.AlignmentFile(src, "rb") as fin:
        with pysam.AlignmentFile(dst, "wb", header=fin.header) as fout:
            for chrom, pos in regions:
                for read in fin.fetch(chrom, max(0, pos - window), pos + window):
                    n_in += 1
                    fout.write(read)
                    n_out += 1
    pysam.index(dst)
    return n_out


def main():
    out = os.path.abspath(sys.argv[1] if len(sys.argv) > 1 else "demo_bundle")
    os.makedirs(out, exist_ok=True)
    lines = []
    for label, spec in BUNDLE.items():
        bam_dst = os.path.join(out, f"{label}.bam")
        print(f"slicing {label} ...", flush=True)
        n = slice_bam(spec["bam"], bam_dst, spec["regions"])
        bcf_src = os.path.expanduser(spec["bcf"])
        bcf_dst = os.path.join(out, f"{label}.bcf")
        with open(bcf_src, "rb") as a, open(bcf_dst, "wb") as b:
            b.write(a.read())
        for ext in (".csi", ".tbi"):
            if os.path.exists(bcf_src + ext):
                with open(bcf_src + ext, "rb") as a, open(bcf_dst + ext, "wb") as b:
                    b.write(a.read())
        lines.append((label, n, spec["regions"], spec["why"]))
        print(f"  {label}: {n} reads, "
              f"{os.path.getsize(bam_dst)/1e6:.1f} MB bam + "
              f"{os.path.getsize(bcf_dst)/1e3:.0f} KB bcf", flush=True)

    conf = os.path.join(out, "sv-assistant.conf")
    with open(conf, "w") as f:
        f.write("# Demo bundle. Paths are relative to THIS file, so the bundle\n"
                "# works wherever it is unzipped.\n[paths]\ndata_dir = .\n\n[datasets]\n")
        for label in BUNDLE:
            f.write(f"{label} = {label}.bam\n")
        f.write("\n[candidates]\n")
        for label in BUNDLE:
            f.write(f"{label} = {label}.bcf\n")

    with open(os.path.join(out, "DEMO.md"), "w") as f:
        f.write("# Demo bundle\n\nStart with:\n\n"
                "    SV_CONFIG=<this folder>/sv-assistant.conf python -m stage1_igv_assistant.ui\n\n"
                "## What is here\n\n")
        for label, n, regions, why in lines:
            f.write(f"- **{label}** — {why}\n")
            f.write(f"  - {n} reads kept, in ±{WINDOW:,} bp windows around:\n")
            for c, p in regions:
                f.write(f"    - `{c}:{p:,}`\n")
        clean = measured(os.path.join(out, "DEMO_CLEAN.bam"), "chr20", 200_000)
        repeat = measured(os.path.join(out, "DEMO_REPEAT.bam"), "chr20", 25_800_000)
        f.write("\n## Loci that mean something\n\n"
                "Measured on these slices when the bundle was built, by breakpoint_evidence_summary "
                "(window 500 bp, the layers applicable_layers finds in 20,000 sampled reads).\n\n"
                "| dataset | coordinate | what it shows |\n|---|---|---|\n"
                f"| DEMO_CLEAN | chr20:200000 | {clean} |\n"
                "| DEMO_CLEAN | chr21:14100000 | the partner side of the same junction |\n"
                f"| DEMO_REPEAT | chr20:25800000 | {repeat} |\n"
                "\n## Limits\n\n"
                "Only the windows above contain reads. Any other coordinate returns zero reads and "
                "the layers report `assessable: false`. That is honest, but it is not what the full "
                "BAM would show. The candidate BCFs are complete, so the filter chain and the "
                "two-sample comparison behave exactly as on the full data.\n")
    total = sum(os.path.getsize(os.path.join(out, f)) for f in os.listdir(out))
    print(f"\nbundle: {out}")
    for f in sorted(os.listdir(out)):
        print(f"  {os.path.getsize(os.path.join(out,f))/1e6:8.2f} MB  {f}")
    print(f"  {total/1e6:8.2f} MB  TOTAL")


if __name__ == "__main__":
    main()
