#!/usr/bin/env python3
"""Build the interface's local gene table from a GENCODE or Ensembl GTF (Phase 26).

    python3 scripts/make_gene_table.py GTF[.gz] OUT.tsv.gz

One line per gene: chrom, start, end, strand, gene_name, gene_id, gene_type,
canonical transcript id, and that transcript's exon starts and ends (1-based,
ascending, comma-separated). The canonical transcript is the one the GTF tags
"Ensembl_canonical"; failing that "MANE_Select"; failing that the transcript
with the most exonic bases. The header records the GTF's file name (never its
directory) and sha256, so a figure can be traced to the annotation behind it.

Suggested source (GRCh38, the reference of every file in this project):
  https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/gencode.v46.basic.annotation.gtf.gz
Optional companions read by the interface (see stage1_igv_assistant/tools/gene_table.py):
  https://omim.org/static/omim/data/mim2gene.txt
  https://purl.obolibrary.org/obo/hp/hpoa/genes_to_disease.txt
"""
import gzip
import hashlib
import os
import re
import sys
import time

ATTR = re.compile(r'(\S+) "([^"]*)"')


def opener(p):
    return gzip.open(p, "rt") if p.endswith(".gz") else open(p)


def sha256(p):
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main(gtf, out):
    genes, order = {}, []
    tx = {}                      # transcript id -> {"gene", "tags", "exons"}
    with opener(gtf) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            feat = f[2]
            if feat not in ("gene", "transcript", "exon"):
                continue
            attrs, tags = {}, set()
            for k, v in ATTR.findall(f[8]):
                if k == "tag":
                    tags.add(v)
                else:
                    attrs.setdefault(k, v)
            gid = attrs.get("gene_id")
            if feat == "gene":
                genes[gid] = {"chrom": f[0], "start": int(f[3]), "end": int(f[4]), "strand": f[6],
                              "name": attrs.get("gene_name") or gid,
                              "type": attrs.get("gene_type") or attrs.get("gene_biotype") or "unknown"}
                order.append(gid)
            elif feat == "transcript":
                tx[attrs["transcript_id"]] = {"gene": gid, "tags": tags, "exons": []}
            else:
                t = tx.get(attrs.get("transcript_id"))
                if t is not None:
                    t["exons"].append((int(f[3]), int(f[4])))
    best = {}
    for tid, t in tx.items():
        rank = (2 if "Ensembl_canonical" in t["tags"] else 1 if "MANE_Select" in t["tags"] else 0,
                sum(e - s + 1 for s, e in t["exons"]))
        cur = best.get(t["gene"])
        if cur is None or rank > cur[0]:
            best[t["gene"]] = (rank, tid)
    with gzip.open(out, "wt") as w:
        w.write(f"# local gene table for the breakpoint interface, built {time.strftime('%Y-%m-%d')}\n")
        w.write(f"# source GTF: {os.path.basename(gtf)} sha256 {sha256(gtf)}\n")
        w.write("# columns: chrom start end strand gene_name gene_id gene_type transcript_id exon_starts exon_ends\n")
        n = 0
        for gid in order:
            g = genes[gid]
            tid = best.get(gid, (None, None))[1]
            ex = sorted(tx[tid]["exons"]) if tid else []
            w.write("\t".join([g["chrom"], str(g["start"]), str(g["end"]), g["strand"], g["name"], gid,
                               g["type"], tid or "", ",".join(str(s) for s, _ in ex),
                               ",".join(str(e) for _, e in ex)]) + "\n")
            n += 1
    print(f"{n} genes written; {sum(1 for g in genes.values() if g['type'] == 'protein_coding')} protein-coding")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2])
