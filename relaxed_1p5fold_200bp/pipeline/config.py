"""Configuration for the 2020 MCM3 publication workflow."""

from __future__ import annotations

from pathlib import Path


RUN_ROOT = Path(__file__).resolve().parents[1]
PROJECT_ROOT = RUN_ROOT.parent
PIPELINE_ROOT = RUN_ROOT / "pipeline"

REFERENCE_DIR = RUN_ROOT / "reference"
RNA_ROOT = RUN_ROOT / "deg_work"
CUTRUN_ROOT = RUN_ROOT / "cutrun_work"
REGULATORY_ROOT = RUN_ROOT / "regulatory_work"

RAW_RNA_DIR = PROJECT_ROOT / "bulkRNAseq_data"
RAW_CUTRUN_DIR = PROJECT_ROOT / "CUT&RUN_data"

GENCODE_GTF_NAME = "gencode.vM25.annotation.gtf"
GENCODE_TRANSCRIPTS_NAME = "gencode.vM25.transcripts.fa"
MM10_FASTA_NAME = "GRCm38.primary_assembly.gencodeM25_contigs.fa"
SACCER3_FASTA_NAME = "sacCer3.fa"

FACTORS = ("MCM3", "NONO", "PSPC1")
COLORS = {"MCM3": "#C4161C", "NONO": "#F39C12", "PSPC1": "#4C78A8", "IgG": "#737B8C"}

RNA_FDR_MAX = 0.05
RNA_LOGFC_MIN = 0.28
RNA_MCM3_MIN_CPM = 1.0
RNA_MCM3_MIN_SAMPLES = 2

MAPQ_MIN = 30
PROMOTER_BP = 1000
PROMOTER_OVERLAP_BP = 200
ATOMIC_SEGMENT_MIN_BP = 250
TOTAL_FACTOR_SUPPORT = 2
N_FACTOR_LIBRARIES = 2
RATIO_MIN = 1.5
PEAK_CANDIDATE_P_MAX = 1e-3
PEAK_P_MAX = 1e-4
PEAK_FE_MIN = 3.0
PEAK_Q_2020_MAX = 0.01

YEAST_SCALE_NUMERATOR = 10_000.0
DISPLAY_MULTIPLIER = 10.0


def final_data(*parts: str) -> Path:
    path = RUN_ROOT.joinpath(*parts)
    path.parent.mkdir(parents=True, exist_ok=True)
    return path
