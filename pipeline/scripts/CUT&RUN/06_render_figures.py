#!/usr/bin/env python3
"""Render and publish the final CUT&RUN figures."""

#Before you run this script, please replace directory placeholders with your own directory

from __future__ import annotations

import argparse
import importlib.util
import os
import shutil
import sys
from pathlib import Path

PIPELINE_ROOT = Path(__file__).resolve().parents[1]
if str(PIPELINE_ROOT) not in sys.path:
    sys.path.insert(0, str(PIPELINE_ROOT))

from config import CUTRUN_ROOT, FACTORS, REFERENCE_DIR, RUN_ROOT


TAG = "q5e2_fe3_min2of2"
HELPERS = Path(__file__).resolve().parent / "figures_rendering"
FIGURE_DATA = CUTRUN_ROOT / "data" / "figure_inputs"
VISUALS = CUTRUN_ROOT / "visuals"
PUBLICATION_FIGURES = VISUALS / "publication_figures"
TRACK_ROOT = CUTRUN_ROOT / "03_bigwig" / "IGV_representation"
TRACKS = {factor: TRACK_ROOT / f"{factor}_mean.bw" for factor in FACTORS}
IGG = TRACK_ROOT / "IgG_mean.bw"


def load_module(name: str):
    path = HELPERS / name
    spec = importlib.util.spec_from_file_location(path.stem, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot import {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def executable(name: str, fallback: str) -> str:
    path = shutil.which(name)
    if path:
        return path
    if Path(fallback).is_file():
        return fallback
    raise FileNotFoundError(f"Required executable not found: {name}")


def prepare_inputs() -> None:
    FIGURE_DATA.mkdir(parents=True, exist_ok=True)
    counts = CUTRUN_ROOT / "data" / "Venn_Peaks_counts.tsv"
    if not counts.is_file():
        raise FileNotFoundError(counts)
    shutil.copy2(counts, FIGURE_DATA / "Venn_Peaks_counts.tsv")


def enable_text_free_rendering() -> None:
    """Suppress all Matplotlib text and append `_noTexts` during figure export."""
    from matplotlib.figure import Figure
    from matplotlib.text import Text

    original_savefig = Figure.savefig

    def savefig_without_text(self, fname, *args, **kwargs):
        path = Path(fname)
        if path.suffix.lower() == ".png":
            fname = path.with_name(f"{path.stem}_noTexts{path.suffix}")
        return original_savefig(self, fname, *args, **kwargs)

    Figure.savefig = savefig_without_text
    Text.draw = lambda self, renderer: None


def configure_primary():
    module = load_module("render_profiles.py")
    module.PROJECT = HELPERS
    module.RUN = RUN_ROOT
    module.CUTRUN = CUTRUN_ROOT
    module.DATA = FIGURE_DATA
    module.VISUALS = VISUALS
    module.PROMOTER = FIGURE_DATA / "promoter_gene_venn"
    module.TAG = TAG
    module.GENE_BED = CUTRUN_ROOT / "data" / "GeneBodies_M25.bed6"
    module.TRACKS = TRACKS
    module.IGG = IGG
    module.COMPUTE_MATRIX = Path(executable("computeMatrix", "/opt/anaconda3/envs/cutrun_env/bin/computeMatrix"))
    module.RSCRIPT = executable("Rscript", "Rscript")
    return module


def configure_distribution():
    module = load_module("render_peak_distribution.py")
    module.PROJECT = HELPERS
    module.RUN = RUN_ROOT
    module.CUTRUN = CUTRUN_ROOT
    module.DATA = FIGURE_DATA
    module.VISUALS = VISUALS
    module.TEXT_FREE = False
    module.PEAKS = FIGURE_DATA / "batch_stratified_peaks"
    module.STAGE = FIGURE_DATA / "peak_distribution" / "intermediate"
    module.GTF = REFERENCE_DIR / "gencode.vM25.annotation.gtf"
    module.IGG_BED = CUTRUN_ROOT / "data" / "PeakSets" / "Peaks_IgG.bed"
    module.PEAK_LOCI = {
        factor: CUTRUN_ROOT / "data" / "PeakLoci" / f"Venn_Peaks_{factor}.bed"
        for factor in FACTORS
    }
    module.BEDTOOLS = Path(executable("bedtools", "/opt/anaconda3/envs/cutrun_env/bin/bedtools"))
    return module


def configure_additional():
    module = load_module("render_peak_profiles.py")
    module.PROJECT = HELPERS
    module.RUN = RUN_ROOT
    module.CUTRUN = CUTRUN_ROOT
    module.DATA = FIGURE_DATA
    module.VISUALS = VISUALS
    module.TAG = TAG
    module.GTF = REFERENCE_DIR / "gencode.vM25.annotation.gtf"
    module.GENE_BED = CUTRUN_ROOT / "data" / "GeneBodies_M25.bed6"
    module.TRACKS = TRACKS
    module.PEAK_LOCI = CUTRUN_ROOT / "data" / "PeakLoci" / "Venn_Peaks_loci.tsv"
    module.PROMOTERS = CUTRUN_ROOT / "data" / "Promoters_M25_TSSplusminus1kb.bed"
    module.OUT = FIGURE_DATA / "additional_panels"
    module.COMPUTE_MATRIX = Path(executable("computeMatrix", "/opt/anaconda3/envs/cutrun_env/bin/computeMatrix"))
    return module


def main() -> None:
    global VISUALS
    parser = argparse.ArgumentParser()
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--rpkm-only", action="store_true")
    parser.add_argument("--publication-figures", action="store_true", help="Also render text-free PNGs in cutrun_work/visuals/publication_figures")
    parser.add_argument("--no-text", action="store_true", help=argparse.SUPPRESS)
    args = parser.parse_args()
    if args.no_text:
        VISUALS = PUBLICATION_FIGURES
        enable_text_free_rendering()
    if args.dry_run:
        print("[DRY-RUN] Would render the registered CUT&RUN figures.")
        return
    prepare_inputs()
    if args.publication_figures:
        PUBLICATION_FIGURES.mkdir(parents=True, exist_ok=True)
    required = [REFERENCE_DIR / "gencode.vM25.annotation.gtf", CUTRUN_ROOT / "data" / "GeneBodies_M25.bed6", IGG, *TRACKS.values()]
    missing = [str(path) for path in required if not path.is_file() or path.stat().st_size == 0]
    if missing:
        raise FileNotFoundError("Missing CUT&RUN figure inputs:\n" + "\n".join(missing))
    primary = configure_primary()
    primary.NO_TEXT_VISUALS = None
    primary.TEXT_FREE = args.no_text
    if args.rpkm_only:
        primary.render_rpkm_profiles(primary.promoter_sets())
        if args.publication_figures and not args.no_text:
            subprocess.run([sys.executable, str(Path(__file__).resolve()), "--no-text", "--rpkm-only"], check=True)
        print(f"[DONE] CUT&RUN RPKM figures: {VISUALS}")
        return
    primary.main()
    distribution = configure_distribution()
    distribution.TEXT_FREE = args.no_text
    distribution.main()
    counts = FIGURE_DATA / "peak_distribution" / "Pie_PeakDistribution_counts.tsv"
    shutil.copy2(counts, CUTRUN_ROOT / "data" / "Pie_PeakDistribution_counts.tsv")
    additional = configure_additional()
    additional.NO_TEXT_VISUALS = None
    additional.main()
    peak_gene_module = load_module(HELPERS / "render_peak_associated_genes.py")
    peak_gene_module.VISUALS = VISUALS
    peak_gene_module.NO_TEXT_VISUALS = None
    peak_gene_module.main()
    if args.publication_figures and not args.no_text:
        subprocess.run([sys.executable, str(Path(__file__).resolve()), "--no-text"], check=True)
        print(f"[DONE] CUT&RUN publication figures: {PUBLICATION_FIGURES}")
    print(f"[DONE] CUT&RUN figures: {VISUALS}")


if __name__ == "__main__":
    for variable in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"):
        os.environ[variable] = "1"
    os.environ.setdefault("MPLCONFIGDIR", str(CUTRUN_ROOT / "data" / ".matplotlib"))
    Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)
    main()
