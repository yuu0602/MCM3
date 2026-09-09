#!/usr/bin/env python3
"""Render and publish the final CUT&RUN figures."""

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


TAG = "p1e4_q5e2_fe3_union4of6"
HELPERS = Path(__file__).resolve().parent / "figures_rendering"
FIGURE_DATA = CUTRUN_ROOT / "data" / "figure_inputs"
VISUALS = CUTRUN_ROOT / "visuals"
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
    module.RSCRIPT = executable("Rscript", "/usr/local/bin/Rscript")
    return module


def configure_distribution():
    module = load_module("render_peak_distribution.py")
    module.PROJECT = HELPERS
    module.RUN = RUN_ROOT
    module.CUTRUN = CUTRUN_ROOT
    module.DATA = FIGURE_DATA
    module.VISUALS = VISUALS
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


def publish(mapping: dict[str, str]) -> None:
    VISUALS.mkdir(parents=True, exist_ok=True)
    for source_name, target_name in mapping.items():
        source, target = VISUALS / source_name, VISUALS / target_name
        if source == target:
            if not target.is_file():
                raise FileNotFoundError(target)
            continue
        if not source.is_file():
            if target.is_file():
                continue
            raise FileNotFoundError(source)
        target.unlink(missing_ok=True)
        source.replace(target)


OUTPUTS = {
    "FigX_TSS_to_TES_metaprofile__NT_KD_UNION__minus3_to_plus3_TES2p8.png": "Metaprofile_PromoterGenes.png",
    "FigS2e_peak_overlap_whole_genome.png": "Venn_Peaks.png",
    "FigS2e_peak_overlap_whole_genome_no_numbers.png": "Venn_Peaks_noNumbers.png",
    "MCM3_RPKM.png": "MCM3_RPKM.png",
    "NONO_RPKM.png": "NONO_RPKM.png",
    "PSPC1_RPKM.png": "PSPC1_RPKM.png",
    "FigX_promoter_only_profile__region_101.png": "Profile_MCM3_PSPC1.png",
    "FigX_promoter_only_profile__region_110.png": "Profile_MCM3_NONO.png",
    "FigX_promoter_only_profile__region_111.png": "Profile_MCM3_NONO_PSPC1.png",
    f"FigX_promoter_gene_overlap_venn__{TAG}.png": "Venn_PromoterGenes.png",
    f"FigX_promoter_gene_overlap_venn__{TAG}_no_numbers.png": "Venn_PromoterGenes_noNumbers.png",
    "FigX_peak_distribution_IgG.png": "Pie_IgG.png",
    "FigX_peak_distribution_MCM3.png": "Pie_MCM3.png",
    "FigX_peak_distribution_NONO.png": "Pie_NONO.png",
    "FigX_peak_distribution_PSPC1.png": "Pie_PSPC1.png",
    "FigX_peak_distribution_pies.png": "Pie_PeakDistribution.png",
    "MCM3_Peak_profiles.png": "MCM3_Peak_profiles.png",
    "MCM3_NONO_PSPC1_peak_profiles.png": "MCM3_NONO_PSPC1_peak_profiles.png",
}

def clean() -> None:
    for path in VISUALS.glob("._*"):
        path.unlink()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--publish-only", action="store_true")
    args = parser.parse_args()
    if args.dry_run:
        print("[DRY-RUN] Would render the registered CUT&RUN figures.")
        return
    if not args.publish_only:
        prepare_inputs()
        required = [REFERENCE_DIR / "gencode.vM25.annotation.gtf", CUTRUN_ROOT / "data" / "GeneBodies_M25.bed6", IGG, *TRACKS.values()]
        missing = [str(path) for path in required if not path.is_file() or path.stat().st_size == 0]
        if missing:
            raise FileNotFoundError("Missing CUT&RUN figure inputs:\n" + "\n".join(missing))
        configure_primary().main()
        distribution = configure_distribution()
        distribution.main()
        counts = FIGURE_DATA / "peak_distribution" / "Pie_PeakDistribution_counts.tsv"
        shutil.copy2(counts, CUTRUN_ROOT / "data" / "Pie_PeakDistribution_counts.tsv")
        configure_additional().main()
    publish(OUTPUTS)
    clean()
    print(f"[DONE] CUT&RUN figures: {VISUALS}")


if __name__ == "__main__":
    for variable in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"):
        os.environ[variable] = "1"
    os.environ.setdefault("MPLCONFIGDIR", str(CUTRUN_ROOT / "data" / ".matplotlib"))
    Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)
    main()
