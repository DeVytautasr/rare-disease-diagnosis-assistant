"""
review_fixture.py -- synthetic reads, calls and genes for the Phase 26 review tests.

Not a test itself (no test_ prefix): test_junction_review.py and the interface
tests build it in a temporary directory. Everything is invented; the chromosome
lengths are GRCh38's, so bands and centromeres resolve as on real data.

Planted rearrangements (every read below is placed on purpose; the counts the
tests expect are the counts planted here):

  R1  balanced t(19;22), both breakpoints on long arms: both junctions called
        J1 3to5  chr19:48,000,000 (left part)  + chr22:30,000,000 (right part)
                 J1_PAIRS read pairs, J1_SPLIT split reads, and J1_ONE_SIDED pair whose
                 mate maps ambiguously (not counted)
        J2 5to3  chr22:29,999,999 (left part)  + chr19:48,000,201 (right part)
                 J2_PAIRS read pairs, J2_SPLIT split reads
        200 bp of chr19 (48,000,001-48,000,200) lie between the two junctions.
        GENE_A (+) is broken in intron 2 of 3; GENE_B (-) in intron 2 of 3.
  R2  one junction chr4:60,000,000 (long arm) <-> chr19:10,000,000 (short arm), 5to5
        (R2_PAIRS, R2_SPLIT): chr4's distal long arm joined to chr19's centric part;
        its reciprocal junction is in the file but LowQual (PE 2, SR 0), and
        R2_RECIP_PAIRS pairs support that reciprocal join in the reads.
  R3  a segment of chr22 (40,000,000-40,050,000) joined into chr4 at 100,000,000:
        two junctions, near on chr4 and 50 kb apart on chr22. GENE_D (OMIM) inside.
  R4  chr4:50,000,000 <-> chr19:26,000,000, both in centromeres: its joining pairs
        are all at mapping quality 0-10, and 20 pairs at the chr4 end point to
        chr1, chr7 and chr16.
  R5  chr19:45,000,000 <-> chr22:20,000,000, called but with no read behind it.
        Also called in the second sample (OTHER), so it is recurrent.
"""
import os
import random

import pysam

LEN = {"chr1": 248956422, "chr4": 190214555, "chr7": 159345973, "chr16": 90338345,
       "chr19": 58617616, "chr22": 50818468}
TID = {c: i for i, c in enumerate(LEN)}
J1_PAIRS, J1_SPLIT = 12, 6
J1_ONE_SIDED = 1          # a J1 pair with its mate at mapping quality 5
J2_PAIRS, J2_SPLIT = 9, 4
R2_PAIRS, R2_SPLIT, R2_RECIP_PAIRS = 10, 5, 4
R3A_PAIRS, R3B_PAIRS = 6, 5
R4_ELSEWHERE = 20
R4_AMBIGUOUS = 80

_rng = random.Random(26)


def _seq(n):
    return "".join(_rng.choice("ACGT") for _ in range(n))


def _seg(h, name, chrom, start1, cigar, flag, mchrom, mstart1, tlen=0, mapq=60, tags=None):
    """start1 / mstart1 are 1-based; cigar is a list of (op, length)."""
    r = pysam.AlignedSegment(h)
    qlen = sum(n for op, n in cigar if op in (0, 1, 4))
    r.query_name = name
    r.query_sequence = _seq(qlen)
    r.query_qualities = pysam.qualitystring_to_array("I" * qlen)
    r.reference_id, r.reference_start, r.mapping_quality = TID[chrom], start1 - 1, mapq
    r.flag, r.cigar = flag, cigar
    if mchrom is None:
        r.next_reference_id, r.next_reference_start = TID[chrom], start1 - 1
    else:
        r.next_reference_id, r.next_reference_start = TID[mchrom], mstart1 - 1
    r.template_length = tlen
    r.set_tag("RG", "rg1")
    for k, v in (tags or {}).items():
        r.set_tag(k, v)
    return r


P, PROPER, R1F, R2F, REV, MREV, SUPP, MUNMAP = 0x1, 0x2, 0x40, 0x80, 0x10, 0x20, 0x800, 0x8


def _pair(h, out, name, c1, s1, rev1, c2, s2, rev2, mapq=60, mq2=None):
    """A discordant pair: read 1 at c1:s1, read 2 at c2:s2 (both 150 M)."""
    f1 = P | R1F | (REV if rev1 else 0) | (MREV if rev2 else 0)
    f2 = P | R2F | (REV if rev2 else 0) | (MREV if rev1 else 0)
    out.append(_seg(h, name, c1, s1, [(0, 150)], f1, c2, s2, 0, mapq, {"MQ": mq2 if mq2 is not None else mapq}))
    out.append(_seg(h, name, c2, s2, [(0, 150)], f2, c1, s1, 0, mq2 if mq2 is not None else mapq, {"MQ": mapq}))


def _split_left_right(h, out, name, cl, end_l, cr, start_r, k):
    """A split read whose first k bases end chromosome cl at end_l (left part) and
    whose remaining 150-k bases start chromosome cr at start_r (right part).
    Primary on cl (k M, then clipped), supplementary on cr; its mate a normal read on cl."""
    s_l = end_l - k + 1
    mate = s_l - 300
    out.append(_seg(h, name, cl, s_l, [(0, k), (4, 150 - k)], P | PROPER | R2F | REV, cl, mate, -(s_l + k - mate), 60,
                    {"SA": f"{cr},{start_r},-,{k}S{150 - k}M,60,0;", "MQ": 60}))
    out.append(_seg(h, name, cl, mate, [(0, 150)], P | PROPER | R1F | MREV, cl, s_l, s_l + k - mate, 60, {"MQ": 60}))
    out.append(_seg(h, name, cr, start_r, [(5, k), (0, 150 - k)], P | R2F | REV | SUPP, cl, mate, 0, 60,
                    {"SA": f"{cl},{s_l},-,{k}M{150 - k}S,60,0;"}))


def _split_right_right(h, out, name, ca, start_a, cb, start_b, k):
    """A split read for a 5to5 join: the first k bases are the reverse of cb's right part,
    the rest cb... modelled as primary on ca starting at start_a (left clipped) and a
    supplementary piece on cb starting at start_b, on the opposite strand, left clipped."""
    mate = start_a + 400
    out.append(_seg(h, name, ca, start_a, [(4, k), (0, 150 - k)], P | PROPER | R1F | MREV, ca, mate, 550, 60,
                    {"SA": f"{cb},{start_b},-,{k}S{150 - k}M,60,0;", "MQ": 60}))
    out.append(_seg(h, name, ca, mate, [(0, 150)], P | PROPER | R2F | REV, ca, start_a, -550, 60, {"MQ": 60}))
    out.append(_seg(h, name, cb, start_b, [(5, k), (0, 150 - k)], P | R1F | REV | SUPP, ca, mate, 0, 60,
                    {"SA": f"{ca},{start_a},+,{k}S{150 - k}M,60,0;"}))


def _background(h, out, chrom, centre, half=3000, depth=30, mapq=60):
    a, b = max(1, centre - half), centre + half
    n = int((b - a) * depth / 300)
    for i in range(n):
        s = _rng.randint(a, b - 560)
        ins = int(_rng.gauss(400, 30))
        m = s + ins - 150
        name = f"bg_{chrom}_{centre}_{i}"
        out.append(_seg(h, name, chrom, s, [(0, 150)], P | PROPER | R1F | MREV, chrom, m, ins, mapq, {"MQ": mapq}))
        out.append(_seg(h, name, chrom, m, [(0, 150)], P | PROPER | R2F | REV, chrom, s, -ins, mapq, {"MQ": mapq}))


def build_reads(path, sample="SYNTH"):
    h = pysam.AlignmentHeader.from_dict({
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": c, "LN": n} for c, n in LEN.items()],
        "RG": [{"ID": "rg1", "SM": sample}]})
    out = []
    for c, p in (("chr19", 48000100), ("chr22", 30000000), ("chr4", 60000000), ("chr19", 10000000),
                 ("chr4", 100000000), ("chr22", 40000000), ("chr22", 40050000), ("chr4", 50000000),
                 ("chr19", 26000000), ("chr19", 45000000), ("chr22", 20000000)):
        _background(h, out, c, p)
    # R1 / J1 3to5: chr19 left part (+ reads) to chr22 right part (- mates)
    for i in range(J1_PAIRS):
        _pair(h, out, f"j1p{i}", "chr19", 48000000 - 460 + 25 * i, False, "chr22", 30000010 + 20 * i, True)
    for i in range(J1_SPLIT):
        _split_left_right(h, out, f"j1s{i}", "chr19", 48000000, "chr22", 30000000, 60 + 8 * i)
    # one more J1 pair whose mate maps ambiguously: seen from chr19 only, so not counted
    _pair(h, out, "j1lowmate", "chr19", 47999860, False, "chr22", 30000420, True, mapq=60, mq2=5)
    # R1 / J2 5to3: chr22 left part (+ reads) to chr19 right part (- mates)
    for i in range(J2_PAIRS):
        _pair(h, out, f"j2p{i}", "chr22", 29999999 - 460 + 30 * i, False, "chr19", 48000210 + 25 * i, True)
    for i in range(J2_SPLIT):
        _split_left_right(h, out, f"j2s{i}", "chr22", 29999999, "chr19", 48000201, 70 + 9 * i)
    # R2 5to5: chr4 right part (- reads) with chr19 right part (- mates)
    for i in range(R2_PAIRS):
        _pair(h, out, f"r2p{i}", "chr4", 60000010 + 22 * i, True, "chr19", 10000015 + 21 * i, True)
    for i in range(R2_SPLIT):
        _split_right_right(h, out, f"r2s{i}", "chr4", 60000000, "chr19", 10000000, 55 + 7 * i)
    # ... and reads for its uncalled reciprocal join, 3to3: + reads at both left parts
    for i in range(R2_RECIP_PAIRS):
        _pair(h, out, f"r2r{i}", "chr4", 60000150 - 470 + 30 * i, False, "chr19", 10000050 - 460 + 35 * i, False)
    # R3: chr4 left + chr22 segment start (3to5), chr22 segment end + chr4 right (5to3)
    for i in range(R3A_PAIRS):
        _pair(h, out, f"r3a{i}", "chr4", 100000000 - 450 + 30 * i, False, "chr22", 40000020 + 25 * i, True)
    for i in range(R3B_PAIRS):
        _pair(h, out, f"r3b{i}", "chr22", 40050000 - 450 + 30 * i, False, "chr4", 100000020 + 25 * i, True)
    _split_left_right(h, out, "r3s0", "chr4", 100000000, "chr22", 40000000, 75)
    # R4: joining pairs only at low mapping quality; other pairs to three other chromosomes
    for i in range(4):
        _pair(h, out, f"r4p{i}", "chr4", 50000000 - 400 + 50 * i, False, "chr19", 26000020 + 30 * i, True,
              mapq=3 + i, mq2=0)
    for i in range(R4_ELSEWHERE):
        oc = ("chr1", "chr7", "chr16")[i % 3]
        _pair(h, out, f"r4e{i}", "chr4", 49999300 + 60 * i, i % 2 == 1, oc, 1000000 + 7919 * i, False, mapq=35)
    for i in range(R4_AMBIGUOUS):            # ambiguous normal reads at the chr19 centromere end
        s = 25999100 + 22 * i
        out.append(_seg(h, f"r4z{i}", "chr19", s, [(0, 150)], P | PROPER | R1F | MREV, "chr19", s + 250, 400, 0))
    out.sort(key=lambda r: (r.reference_id, r.reference_start))
    with pysam.AlignmentFile(path, "wb", header=h) as f:
        for r in out:
            f.write(r)
    pysam.index(path)
    return path


VCF_HEAD = "".join([
    "##fileformat=VCFv4.2\n",
    *[f"##contig=<ID={c},length={n}>\n" for c, n in LEN.items()],
    '##FILTER=<ID=LowQual,Description="Poor quality">\n',
    '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="type">\n',
    '##INFO=<ID=CHR2,Number=1,Type=String,Description="chr2">\n',
    '##INFO=<ID=POS2,Number=1,Type=Integer,Description="pos2">\n',
    '##INFO=<ID=END,Number=1,Type=Integer,Description="end">\n',
    '##INFO=<ID=CT,Number=1,Type=String,Description="connection">\n',
    '##INFO=<ID=PE,Number=1,Type=Integer,Description="pairs">\n',
    '##INFO=<ID=SR,Number=1,Type=Integer,Description="split">\n',
    '##INFO=<ID=PRECISE,Number=0,Type=Flag,Description="precise">\n',
    '##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="imprecise">\n',
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"])


def _bnd(c, p, rid, c2, p2, ct, pe, sr, filt="PASS"):
    return (f"{c}\t{p}\t{rid}\tN\t<TRA>\t.\t{filt}\t"
            f"SVTYPE=BND;CHR2={c2};POS2={p2};CT={ct};PE={pe};SR={sr};PRECISE\n")


def build_calls(path, other=False):
    """The caller's file. other=True writes the second sample's: R5 and noise only."""
    recs = []
    if not other:
        recs += [_bnd("chr19", 48000000, "BND00000001", "chr22", 30000000, "3to5", 14, 6),
                 _bnd("chr19", 48000201, "BND00000002", "chr22", 29999999, "5to3", 10, 4),
                 _bnd("chr4", 60000000, "BND00000003", "chr19", 10000000, "5to5", 11, 5),
                 _bnd("chr4", 60000150, "BND00000004", "chr19", 10000050, "3to3", 2, 0, "LowQual"),
                 _bnd("chr4", 100000000, "BND00000005", "chr22", 40000000, "3to5", 7, 1),
                 _bnd("chr4", 100000001, "BND00000006", "chr22", 40050000, "5to3", 6, 1),
                 _bnd("chr4", 50000000, "BND00000007", "chr19", 26000000, "3to5", 4, 1),
                 "chr19\t50000000\tDEL00000008\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=50010000;CT=3to5;PE=6;SR=2;PRECISE\n"]
    recs.append(_bnd("chr19", 45000000, "BND00000009", "chr22", 20000000, "3to5", 3, 1))
    for i in range(12):                       # noise the filters remove
        recs.append(_bnd("chr16", 1000000 + 50000 * i, f"BND0000010{i:02d}", "chr7", 2000000 + 40000 * i,
                         "3to5", i % 3, 0, "LowQual" if i % 2 else "PASS"))
    with open(path, "w") as f:
        f.write(VCF_HEAD)
        for r in sorted(recs, key=lambda x: (list(LEN).index(x.split("\t")[0]), int(x.split("\t")[1]))):
            f.write(r)
    return path


GENES = [
    # chrom, start, end, strand, name, id, type, transcript, exon starts, exon ends
    ("chr19", 47990000, 48050000, "+", "GENE_A", "ENSG_A", "protein_coding", "ENST_A",
     "47990000,47995000,48010000,48049000", "47990500,47995200,48010300,48050000"),
    ("chr19", 48000050, 48000600, "-", "LNC_F", "ENSG_F", "lncRNA", "ENST_F", "48000050", "48000600"),
    ("chr22", 29950000, 30100000, "-", "GENE_B", "ENSG_B", "protein_coding", "ENST_B",
     "29950000,29990000,30020000,30099000", "29951000,29990200,30020100,30100000"),
    ("chr4", 59950000, 60100000, "+", "GENE_C", "ENSG_C", "protein_coding", "ENST_C",
     "59950000,59980000,60050000,60099500", "59950400,59980150,60050300,60100000"),
    ("chr22", 40010000, 40030000, "-", "GENE_D", "ENSG_D", "protein_coding", "ENST_D",
     "40010000,40029000", "40011000,40030000"),
    ("chr19", 10100000, 10200000, "+", "GENE_E", "ENSG_E", "protein_coding", "ENST_E",
     "10100000,10199000", "10100500,10200000"),
    ("chr22", 45000000, 45100000, "+", "GENE_G", "ENSG_G", "protein_coding", "ENST_G",
     "45000000", "45100000"),
]


def build_genes(directory):
    g = os.path.join(directory, "genes.tsv")
    with open(g, "w") as f:
        f.write("# synthetic gene table for the Phase 26 tests\n")
        for row in GENES:
            f.write("\t".join(str(x) for x in row) + "\n")
    m = os.path.join(directory, "mim2gene.txt")
    with open(m, "w") as f:
        f.write("# MIM Number\tMIM Entry Type\tEntrez Gene ID\tApproved Gene Symbol\tEnsembl Gene ID\n")
        f.write("600001\tgene/phenotype\t1\tGENE_A\tENSG_A\n600004\tgene\t4\tGENE_D\tENSG_D\n"
                "600005\tgene\t5\tGENE_E\tENSG_E\n")
    d = os.path.join(directory, "genes_to_disease.txt")
    with open(d, "w") as f:
        f.write("ncbi_gene_id\tgene_symbol\tassociation_type\tdisease_id\tsource\n")
        f.write("NCBIGene:1\tGENE_A\tMENDELIAN\tOMIM:610002\tx\nNCBIGene:1\tGENE_A\tMENDELIAN\tORPHA:123\tx\n"
                "NCBIGene:4\tGENE_D\tMENDELIAN\tOMIM:610001\tx\n")
    return g, m, d


def build_mask(directory):
    p = os.path.join(directory, "excl.tsv")
    with open(p, "w") as f:
        f.write("chr16\t1000000\t1100000\tsynthetic\n")
    return p


def build_all(directory):
    os.makedirs(directory, exist_ok=True)
    bam = build_reads(os.path.join(directory, "synth.bam"))
    calls = build_calls(os.path.join(directory, "synth.vcf"))
    other = build_calls(os.path.join(directory, "other.vcf"), other=True)
    genes = build_genes(directory)
    mask = build_mask(directory)
    return {"bam": bam, "calls": calls, "other": other, "genes": genes, "mask": mask}
