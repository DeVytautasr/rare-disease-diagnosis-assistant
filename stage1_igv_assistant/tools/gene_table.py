"""
gene_table.py -- genes from a table on this computer (Phase 26).

The page and the assistant name the genes a rearrangement breaks, or carries
with the piece that moved, without sending a position anywhere: the lookup
reads a local file. (gene_at_locus, in server.py, asks Ensembl over the
internet; the interface no longer offers it once this table is set up.)

The table is built once from a GENCODE or Ensembl GTF by
scripts/make_gene_table.py. Two optional companions mark the genes OMIM links
to a disorder, which is what a clinical geneticist looks for first:

  mim2gene.txt          OMIM's free gene list: the gene's MIM number, and
                        whether the entry also carries a phenotype
                        ("gene/phenotype")
  genes_to_disease.txt  the HPO project's gene -> disorder table (OMIM and
                        Orphanet identifiers)

Table format (tab-separated, lines starting with '#' are comments):
  chrom start end strand gene_name gene_id gene_type transcript_id exon_starts exon_ends
1-based inclusive coordinates; the exons are those of the gene's canonical
transcript, ascending, comma-separated (empty when the GTF named none).
"""

import bisect as _bisect
import gzip as _gzip
import os as _os

from stage1_igv_assistant.tools.junction_tools import norm_chrom, display_chrom

NEAREST_MAX_BP = 1_000_000      # author judgement: nearest genes are looked for within 1 Mb
SEGMENT_LIST_LIMIT = 200        # genes listed for one segment; the count is never capped


def _open(path):
    return _gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


class GeneTable:
    def __init__(self, path, mim2gene_path=None, disease_path=None):
        self.path, self.error, self.header = path, None, []
        self.by_chrom = {}            # norm chrom -> list of genes sorted by start
        self.starts = {}
        self.max_len = {}
        self.n_genes = 0
        self.mim = {}                 # symbol -> (gene MIM, has phenotype)
        self.disorders = {}           # symbol -> sorted disorder ids
        self.sources = {}
        try:
            with _open(path) as fh:
                for line in fh:
                    if line.startswith("#"):
                        self.header.append(line[1:].strip())
                        continue
                    f = line.rstrip("\n").split("\t")
                    if len(f) < 7:
                        continue
                    g = {"chromosome": display_chrom(f[0]), "start": int(f[1]), "end": int(f[2]),
                         "strand": f[3], "name": f[4], "gene_id": f[5], "type": f[6],
                         "transcript_id": f[7] if len(f) > 7 and f[7] else None,
                         "exons": _exons(f[8] if len(f) > 8 else "", f[9] if len(f) > 9 else "")}
                    self.by_chrom.setdefault(norm_chrom(f[0]), []).append(g)
                    self.n_genes += 1
        except (OSError, ValueError, IndexError) as e:
            self.error = f"the gene table could not be read ({type(e).__name__})"
            self.by_chrom = {}
            self.n_genes = 0
            return
        for c, genes in self.by_chrom.items():
            genes.sort(key=lambda g: (g["start"], g["end"]))
            self.starts[c] = [g["start"] for g in genes]
            self.max_len[c] = max(g["end"] - g["start"] + 1 for g in genes)
        self.sources["genes"] = self.header[:3]
        if mim2gene_path:
            self._load_mim2gene(mim2gene_path)
        if disease_path:
            self._load_disorders(disease_path)

    # ── optional disease marks ───────────────────────────────────────────────
    def _load_mim2gene(self, path):
        try:
            with _open(path) as fh:
                for line in fh:
                    if line.startswith("#"):
                        continue
                    f = line.rstrip("\n").split("\t")
                    if len(f) < 4 or not f[3]:
                        continue
                    if f[1] in ("gene", "gene/phenotype"):
                        self.mim[f[3]] = (int(f[0]), f[1] == "gene/phenotype")
            self.sources["mim2gene"] = f"{len(self.mim)} genes with an OMIM gene entry"
        except (OSError, ValueError) as e:
            self.sources["mim2gene"] = f"could not be read ({type(e).__name__})"

    def _load_disorders(self, path):
        try:
            d = {}
            with _open(path) as fh:
                for line in fh:
                    f = line.rstrip("\n").split("\t")
                    if len(f) < 4 or f[0].startswith(("ncbi_gene_id", "#")):
                        continue
                    sym, dis = f[1], f[3]
                    if dis.startswith(("OMIM:", "ORPHA:")):
                        d.setdefault(sym, set()).add(dis)
            self.disorders = {k: sorted(v) for k, v in d.items()}
            self.sources["disorders"] = f"{len(self.disorders)} genes linked to a disorder"
        except (OSError, ValueError) as e:
            self.sources["disorders"] = f"could not be read ({type(e).__name__})"

    def marks(self, symbol):
        """OMIM and disorder marks for a gene symbol (empty when none is known)."""
        out = {}
        m = self.mim.get(symbol)
        if m:
            out["omim_gene"] = m[0]
            out["omim_phenotype"] = m[1]
        dis = self.disorders.get(symbol)
        if dis:
            out["disorders"] = dis[:12]
            out["disorder_count"] = len(dis)
            out["omim_disorders"] = sum(1 for x in dis if x.startswith("OMIM:"))
        return out

    @property
    def has_disease_marks(self):
        return bool(self.mim or self.disorders)

    # ── lookups ──────────────────────────────────────────────────────────────
    def _candidates(self, chrom, start, end):
        c = norm_chrom(chrom)
        genes, starts = self.by_chrom.get(c), self.starts.get(c)
        if not genes:
            return []
        lo = _bisect.bisect_left(starts, start - self.max_len[c])
        hi = _bisect.bisect_right(starts, end)
        return [g for g in genes[lo:hi] if g["end"] >= start and g["start"] <= end]

    def _public(self, g, where=None):
        out = {k: g[k] for k in ("name", "gene_id", "type", "strand", "start", "end", "chromosome")}
        if where:
            out["where"] = where
        out.update(self.marks(g["name"]))
        return out

    def at(self, chrom, pos):
        """Genes a breakpoint falls in, each with where in its canonical transcript."""
        found = self._candidates(chrom, pos, pos)
        found.sort(key=lambda g: (g["type"] != "protein_coding", g["start"]))
        return [self._public(g, _where(g, pos)) for g in found]

    def overlapping(self, chrom, start, end, limit=SEGMENT_LIST_LIMIT):
        """Genes overlapping a segment, protein-coding first; with the total count."""
        found = self._candidates(chrom, start, end)
        found.sort(key=lambda g: (g["type"] != "protein_coding", g["start"]))
        coding = sum(1 for g in found if g["type"] == "protein_coding")
        marked = [g for g in found if self.marks(g["name"]).get("disorders") or
                  self.marks(g["name"]).get("omim_phenotype")]
        return {"count": len(found), "protein_coding": coding,
                "with_disorder": len(marked),
                "genes": [self._public(g) for g in found[:limit]],
                "truncated": len(found) > limit}

    def track(self, chrom, start, end, limit=40):
        """Genes overlapping a short window, with their canonical exons, for drawing."""
        found = self._candidates(chrom, start, end)
        found.sort(key=lambda g: (g["type"] != "protein_coding", g["start"]))
        out = []
        for g in found[:limit]:
            d = self._public(g)
            d["exons"] = [[s, e] for s, e in g["exons"] if e >= start and s <= end]
            out.append(d)
        return out

    def nearest(self, chrom, pos, max_distance=NEAREST_MAX_BP):
        """The nearest gene on each side of a position that falls in none."""
        c = norm_chrom(chrom)
        genes = self.by_chrom.get(c) or []
        left = right = None
        for g in self._candidates(chrom, pos - max_distance, pos + max_distance):
            if g["end"] < pos:
                d = pos - g["end"]
                if left is None or d < left[0]:
                    left = (d, g)
            elif g["start"] > pos:
                d = g["start"] - pos
                if right is None or d < right[0]:
                    right = (d, g)
        out = {}
        if left:
            out["left"] = {**self._public(left[1]), "distance_bp": left[0]}
        if right:
            out["right"] = {**self._public(right[1]), "distance_bp": right[0]}
        return out if genes else {}


def _exons(starts, ends):
    if not starts or not ends:
        return []
    try:
        s = [int(x) for x in starts.strip(",").split(",") if x]
        e = [int(x) for x in ends.strip(",").split(",") if x]
    except ValueError:
        return []
    return sorted(zip(s, e))


def _where(g, pos):
    """'exon 3 of 12', 'intron 2 of 11' (numbered in the direction of
    transcription), or 'outside its main transcript'."""
    ex = g["exons"]
    if not ex:
        return "inside the gene"
    n = len(ex)
    if pos < ex[0][0] or pos > ex[-1][1]:
        return "inside the gene, outside its main transcript"
    for i, (s, e) in enumerate(ex):
        if s <= pos <= e:
            k = i + 1 if g["strand"] != "-" else n - i
            return f"exon {k} of {n}"
        if i + 1 < n and e < pos < ex[i + 1][0]:
            k = i + 1 if g["strand"] != "-" else n - i - 1
            return f"intron {k} of {n - 1}"
    return "inside the gene"


_TABLE = None


def load(path, mim2gene_path=None, disease_path=None):
    """Load (or replace) the table this process uses. Returns it, or None."""
    global _TABLE
    if not path or not _os.path.isfile(path):
        _TABLE = None
        return None
    _TABLE = GeneTable(path, mim2gene_path if mim2gene_path and _os.path.isfile(mim2gene_path) else None,
                       disease_path if disease_path and _os.path.isfile(disease_path) else None)
    if _TABLE.error:
        return _TABLE
    return _TABLE


def table():
    return _TABLE if _TABLE is not None and not _TABLE.error else None


def status():
    """What the page's status line and the bootstrap report: never a path."""
    t = _TABLE
    if t is None:
        return {"available": False, "reason": "no gene table is set up"}
    if t.error:
        return {"available": False, "reason": t.error}
    return {"available": True, "genes": t.n_genes, "source": t.sources.get("genes"),
            "omim": t.sources.get("mim2gene"), "disorders": t.sources.get("disorders"),
            "disease_marks": t.has_disease_marks}
