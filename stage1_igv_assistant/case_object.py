"""
case_object.py
Structured case record for the IGV breakpoint assistant.
Inspired by GA4GH Phenopackets. Stores all evidence layers
and the final interpretation in one JSON-serialisable object.
"""

from dataclasses import dataclass, field, asdict
from typing import Optional
from datetime import datetime
import json


@dataclass
class BreakpointEvidence:
    """Evidence collected at one candidate breakpoint."""
    label: str
    chromosome: str
    position: int

    # Evidence layers — None = not yet collected
    locus_stats: Optional[dict] = None
    discordant_pairs: Optional[dict] = None
    soft_clipped_reads: Optional[dict] = None
    split_reads: Optional[dict] = None
    depth_profile: Optional[dict] = None
    gene_at_locus: Optional[dict] = None
    reciprocal_check: Optional[dict] = None

    # Integrated result
    evidence_strength: Optional[str] = None
    evidence_score: Optional[float] = None
    signal_layers: Optional[str] = None
    supporting_observations: list = field(default_factory=list)
    interpretation: Optional[str] = None


@dataclass
class SequencingInfo:
    """Metadata about the BAM file."""
    bam_path: str
    genome_build: str = "GRCh38"
    technology: str = "unknown"  # "short_read_illumina" | "pacbio_hifi" | "ont_nanopore"
    is_paired: bool = True
    mean_coverage: Optional[float] = None
    aligner: Optional[str] = None

    @property
    def has_sa_tags(self) -> bool:
        """Short-read modern alignments and long reads typically have SA tags."""
        return self.aligner in ("bwa-mem", "bwa-mem2", "minimap2") or \
               self.technology in ("pacbio_hifi", "ont_nanopore")

    @property
    def applicable_evidence_layers(self) -> list:
        """Which evidence layers are informative for this technology."""
        layers = ["soft_clipped_reads", "read_depth"]
        if self.is_paired:
            layers.append("discordant_pairs")
        if self.has_sa_tags:
            layers.append("split_reads")
        return layers


@dataclass
class ClinicalInfo:
    """Clinical information about the patient/sample."""
    sample_id: str
    organism: str = "homo_sapiens"
    sex: Optional[str] = None          # "male" | "female" | "unknown"
    age_group: Optional[str] = None    # "pediatric" | "adult"
    hpo_terms: list = field(default_factory=list)      # e.g. ["HP:0001250", "HP:0000486"]
    absent_hpo_terms: list = field(default_factory=list)
    suspected_sv_type: Optional[str] = None  # "balanced_translocation" | "deletion" | "inversion"
    cytogenetic_result: Optional[str] = None # e.g. "t(1;8)(p36;q24)"
    karyotype_chromosomes: Optional[list] = None  # e.g. ["chr1", "chr8"]


@dataclass
class BamCase:
    """
    Top-level case object for one structural variant investigation.

    Usage:
        case = BamCase.new("patient_001", "patient_001.bam")
        case.clinical.hpo_terms = ["HP:0001250"]
        case.clinical.suspected_sv_type = "balanced_translocation"
        case.add_breakpoint("BP1", "chr1", 145507200)
        # ... run tools, populate evidence ...
        case.save("patient_001_case.json")
    """
    case_id: str
    clinical: ClinicalInfo
    sequencing: SequencingInfo
    breakpoints: list = field(default_factory=list)  # list of BreakpointEvidence
    created_at: str = field(default_factory=lambda: datetime.now().isoformat())
    stage: str = "stage1_igv_inspection"
    notes: Optional[str] = None

    @classmethod
    def new(cls, case_id: str, bam_path: str, genome_build: str = "GRCh38") -> "BamCase":
        """Create a new empty case."""
        return cls(
            case_id=case_id,
            clinical=ClinicalInfo(sample_id=case_id),
            sequencing=SequencingInfo(bam_path=bam_path, genome_build=genome_build),
        )

    def add_breakpoint(self, label: str, chromosome: str, position: int) -> BreakpointEvidence:
        """Add a candidate breakpoint to investigate."""
        bp = BreakpointEvidence(label=label, chromosome=chromosome, position=position)
        self.breakpoints.append(bp)
        return bp

    def get_breakpoint(self, label: str) -> Optional[BreakpointEvidence]:
        """Retrieve a breakpoint by label."""
        for bp in self.breakpoints:
            if bp.label == label:
                return bp
        return None

    def to_dict(self) -> dict:
        return asdict(self)

    def save(self, path: str):
        """Save case to JSON."""
        with open(path, "w") as f:
            json.dump(self.to_dict(), f, indent=2)
        print(f"Case saved: {path}")

    @classmethod
    def load(cls, path: str) -> "BamCase":
        """Load case from JSON."""
        with open(path) as f:
            data = json.load(f)
        # Reconstruct nested dataclasses
        data["clinical"] = ClinicalInfo(**data["clinical"])
        data["sequencing"] = SequencingInfo(**data["sequencing"])
        data["breakpoints"] = [BreakpointEvidence(**bp) for bp in data["breakpoints"]]
        return cls(**data)

    def summary(self) -> str:
        """Human-readable summary for reporting."""
        lines = [
            f"Case: {self.case_id}",
            f"BAM: {self.sequencing.bam_path}",
            f"Technology: {self.sequencing.technology}",
            f"Suspected SV: {self.clinical.suspected_sv_type}",
            f"Cytogenetics: {self.clinical.cytogenetic_result}",
            f"Breakpoints investigated: {len(self.breakpoints)}",
        ]
        for bp in self.breakpoints:
            lines.append(
                f"  {bp.label}: {bp.chromosome}:{bp.position} — "
                f"{bp.evidence_strength or 'not yet evaluated'} "
                f"({bp.signal_layers or '?/4'})"
            )
        return "\n".join(lines)


# ── Quick test ───────────────────────────────────────────────────────────────
if __name__ == "__main__":
    case = BamCase.new("test_001", "patient.bam")
    case.clinical.hpo_terms = ["HP:0001250", "HP:0000486"]
    case.clinical.suspected_sv_type = "balanced_translocation"
    case.clinical.cytogenetic_result = "t(1;8)(p36;q24)"
    case.clinical.karyotype_chromosomes = ["chr1", "chr8"]
    case.sequencing.technology = "short_read_illumina"

    bp1 = case.add_breakpoint("BP1", "chr1", 145507200)
    bp1.evidence_strength = "STRONG"
    bp1.evidence_score = 87.5
    bp1.signal_layers = "3/4"
    bp1.supporting_observations = [
        "12 discordant pairs (all mates on chr8)",
        "6 soft-clipped reads (consensus position 145507198)",
        "4 split reads (partner loci on chr8:47000000)",
    ]

    bp2 = case.add_breakpoint("BP2", "chr8", 47000000)
    bp2.evidence_strength = "MODERATE"
    bp2.evidence_score = 62.5
    bp2.signal_layers = "2/4"

    print(case.summary())
    print()

    case.save("/tmp/test_case.json")
    loaded = BamCase.load("/tmp/test_case.json")
    print(f"Loaded case ID: {loaded.case_id}")
    print(f"HPO terms: {loaded.clinical.hpo_terms}")
    print("case_object.py OK")
