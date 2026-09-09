"""Pinned configuration for the final MCM3 publication workflow.

The default execution materializes accepted computational intermediates into
this run directory and renders all result tables and figures locally.  The
optional raw mode rebuilds those intermediates from the original FASTQs.
"""

from __future__ import annotations

from pathlib import Path


RUN_ROOT = Path(__file__).resolve().parents[1]
PROJECT_ROOT = RUN_ROOT.parent
PIPELINE_ROOT = RUN_ROOT / "pipeline"

REFERENCE_DIR = RUN_ROOT / "reference"
RNA_ROOT = RUN_ROOT / "deg_work"
CUTRUN_ROOT = RUN_ROOT / "cutrun_work"
REGULATORY_ROOT = RUN_ROOT / "regulatory_work"

RAW_RNASETS = {
    "20161213": PROJECT_ROOT / "bulkRNAseq_data" / "20161213_Mcm3KD_RNAseq",
    "201709": PROJECT_ROOT / "bulkRNAseq_data" / "201709_Mcm3KD_RNAseq",
    "20180718": PROJECT_ROOT / "bulkRNAseq_data" / "20180718_NonoPspc1KD_RNAseq",
}
RAW_CUTRUN_2019 = PROJECT_ROOT / "CUT&RUN_data" / "20191012"
RAW_CUTRUN_2020 = {
    "20200923": PROJECT_ROOT / "CUT&RUN_data" / "20200923",
    "20200929": PROJECT_ROOT / "CUT&RUN_data" / "20200929",
}

# Accepted intermediate products used by the fast, publication-figure path.
# `01_prepare.py` materializes these once into RUN_ROOT; no later
# stage reads the source locations below.
ACCEPTED_INTERMEDIATE_SOURCES = {
    "reference_gtf": PROJECT_ROOT / "reanalysis2" / "ref" / "gencode.vM25.annotation.gtf",
    "rna_salmon": PROJECT_ROOT / "reanalysis2" / "deg_work" / "salmon",
    "cutrun_2019_bowtie2": PROJECT_ROOT / "reanalysis7(2020IgG)" / "cutrun_work_20191012" / "01_bowtie2",
    "cutrun_2019_bam": PROJECT_ROOT / "reanalysis7(2020IgG)" / "cutrun_work_20191012" / "02_bam_mapq30",
    "cutrun_2019_peaks": PROJECT_ROOT / "reanalysis20" / "two_2020IgG_2019factors_current_thresholds" / "union_support2of4_q0p05" / "cutrun_work_20191012" / "04_macs2_p1e4_q5e2_fe3_min2of4",
    "cutrun_2020_bowtie2": PROJECT_ROOT / "reanalysis11" / "cutrun_work_20200923_20200929" / "01_bowtie2",
    "cutrun_2020_bam": PROJECT_ROOT / "reanalysis11" / "cutrun_work_20200923_20200929" / "02_bam",
    "cutrun_2020_peaks": PROJECT_ROOT / "reanalysis11" / "cutrun_work_20200923_20200929" / "04_macs3_p1e4_q1e2_fe3_min2of2",
    "cutrun_bigwig": PROJECT_ROOT / "reanalysis23" / "six_library_union_support4of6" / "cutrun_work_2019_2020" / "03_bigwig",
    "cutrun_membership": PROJECT_ROOT / "reanalysis23" / "six_library_union_support4of6" / "cutrun_work_2019_2020" / "07_paperfigs" / "data",
}

GENCODE_GTF_NAME = "gencode.vM25.annotation.gtf"
GENCODE_TRANSCRIPTS_NAME = "gencode.vM25.transcripts.fa"
MM10_FASTA_NAME = "GRCm38.primary_assembly.gencodeM25_contigs.fa"
SACCER3_FASTA_NAME = "sacCer3.fa"

FACTORS = ("MCM3", "NONO", "PSPC1")
COLORS = {"MCM3": "#C4161C", "NONO": "#F39C12", "PSPC1": "#4C78A8", "IgG": "#737B8C"}

# RNA-seq analysis criteria.
RNA_FDR_MAX = 0.05
RNA_LOGFC_MIN = 0.28
RNA_MCM3_MIN_CPM = 1.0
RNA_MCM3_MIN_SAMPLES = 2

# Final CUT&RUN membership criteria.
MAPQ_MIN = 30
PROMOTER_BP = 1000
PROMOTER_OVERLAP_BP = 250
TOTAL_FACTOR_SUPPORT = 4
N_FACTOR_LIBRARIES = 6
RATIO_MIN = 2.0
PEAK_P_MAX = 1e-4
PEAK_FE_MIN = 3.0
PEAK_Q_2019_MAX = 0.05
PEAK_Q_2020_MAX = 0.01

# All libraries initially use mouse coverage per 10,000 retained yeast pairs.
YEAST_SCALE_NUMERATOR = 10_000.0

# Display-only batch factors; they do not enter peak membership or gene sets.
DISPLAY_2020_TO_2019 = {"MCM3": 0.556, "NONO": 0.515, "PSPC1": 0.537, "IgG": 0.537355}
DISPLAY_MULTIPLIER = 10.0


def final_data(*parts: str) -> Path:
    path = RUN_ROOT.joinpath(*parts)
    path.parent.mkdir(parents=True, exist_ok=True)
    return path
