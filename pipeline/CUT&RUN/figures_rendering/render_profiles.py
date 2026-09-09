#!/usr/bin/env python3
"""Render promoter-overlap, promoter-profile, and RPKM figures."""

from __future__ import annotations

import importlib.util
import math
from pathlib import Path
import subprocess

import numpy as np
import pandas as pd


PROJECT = Path(__file__).resolve().parent
RUN = PROJECT.parents[2]
CUTRUN = RUN / "cutrun_work"
DATA = CUTRUN / "data" / "figure_inputs"
VISUALS = CUTRUN / "visuals"
PROMOTER = DATA / "promoter_gene_venn"
TAG = "p1e4_q5e2_fe3_union4of6"
FACTORS = ("MCM3", "NONO", "PSPC1")
COLORS = {"MCM3": "#A80F14", "NONO": "#F39C12", "PSPC1": "#1F4AA8", "IgG": "#6B7280"}
GENE_BED = CUTRUN / "data" / "GeneBodies_M25.bed6"
TRACKS = {
    factor: CUTRUN / "03_bigwig" / "IGV_representation" / f"{factor}_mean.bw"
    for factor in FACTORS
}
IGG = CUTRUN / "03_bigwig" / "IGV_representation" / "IgG_mean.bw"
COMPUTE_MATRIX = Path("computeMatrix")
RSCRIPT = "Rscript"


def load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot import {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def read_genes(path: Path) -> set[str]:
    table = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    return set(table["gene"].str.strip()) - {""}


def promoter_sets() -> dict[str, set[str]]:
    return {
        factor: read_genes(PROMOTER / f"Venn_PromoterGenes_{factor}.tsv")
        for factor in FACTORS
    }


def promoter_regions(sets: dict[str, set[str]]) -> dict[str, int]:
    a, b, c = (sets[factor] for factor in FACTORS)
    return {
        "111_MCM3_NONO_PSPC1": len(a & b & c),
        "110_MCM3_NONO": len((a & b) - c),
        "101_MCM3_PSPC1": len((a & c) - b),
        "100_only_MCM3": len(a - b - c),
    }


def render_promoter_venn(sets: dict[str, set[str]]) -> None:
    module = load_module(
        PROJECT / "promoter_venn.py",
        "promoter_venn",
    )
    module.PLOT_TITLE = "CUT&RUN promoter-gene overlap"
    module.plot_venn(
        sets,
        VISUALS / f"FigX_promoter_gene_overlap_venn__{TAG}.png",
        show_numbers=True,
    )
    module.plot_venn(
        sets,
        VISUALS / f"FigX_promoter_gene_overlap_venn__{TAG}_no_numbers.png",
        show_numbers=False,
    )


def render_promoter_profiles(sets: dict[str, set[str]]) -> None:
    module = load_module(
        PROJECT / "promoter_profiles.py",
        "promoter_profiles",
    )
    counts = promoter_regions(sets)
    module.OUTDIR = VISUALS
    module.TMPDIR = DATA / "promoter_profiles" / "intermediate"
    module.MATRIX_TMPDIR = module.TMPDIR / "matrices"
    module.REGION_TABLE = PROMOTER / "Venn_PromoterGenes_regions.tsv"
    module.GENE_BED = GENE_BED
    module.BIGWIGS = TRACKS
    module.COMPUTE_MATRIX = COMPUTE_MATRIX
    module.DISPLAY_SCALE = 1.0
    module.SIGNAL_LABEL = "Coverage"
    module.REGIONS = {
        "111": {
            "data_stem": "Profile_MCM3_NONO_PSPC1",
            "table_region": "111_MCM3_NONO_PSPC1",
            "label": "MCM3-PSPC1-NONO co-bound promoter genes",
            "display_title": "MCM3-PSPC1-NONO Co-localized\ngenes",
            "tracks": FACTORS,
            "expected_n": counts["111_MCM3_NONO_PSPC1"],
        },
        "110": {
            "data_stem": "Profile_MCM3_NONO",
            "table_region": "110_MCM3_NONO",
            "label": "MCM3-NONO-only promoter genes",
            "display_title": "MCM3-NONO-only\npromoter genes",
            "tracks": ("MCM3", "NONO"),
            "expected_n": counts["110_MCM3_NONO"],
        },
        "101": {
            "data_stem": "Profile_MCM3_PSPC1",
            "table_region": "101_MCM3_PSPC1",
            "label": "MCM3-PSPC1-only promoter genes",
            "display_title": "MCM3-PSPC1-only\npromoter genes",
            "tracks": ("MCM3", "PSPC1"),
            "expected_n": counts["101_MCM3_PSPC1"],
        },
    }
    module.main(reuse_matrices=False)
    source = VISUALS / "Profile_PromoterGenes_summary.tsv"
    if source.is_file():
        source.replace(DATA / "promoter_profiles" / "Profile_PromoterGenes_summary.tsv")


def render_rpkm_profiles(sets: dict[str, set[str]]) -> None:
    module = load_module(PROJECT / "rpkm_profiles.py", "rpkm_profiles")
    counts = promoter_regions(sets)
    outdir = DATA / "rpkm_profiles" / "intermediate"
    module.CUTRUN_ROOT = CUTRUN
    module.OUTDIR = outdir
    module.TMPDIR = outdir / "regions"
    module.MATRIX_TMPDIR = outdir / "matrices"
    module.REGION_TABLE = PROMOTER / "Venn_PromoterGenes_regions.tsv"
    module.GENE_BED = GENE_BED
    module.BIGWIGS = TRACKS
    module.COMPUTE_MATRIX = COMPUTE_MATRIX
    module.DISPLAY_SCALE = 1.0
    module.SIGNAL_LABEL = "Log2(RPKM)"
    module.OUTPUT = outdir / "MCM3_RPKM.png"
    module.OUTPUT_NONO = outdir / "NONO_RPKM.png"
    module.OUTPUT_PSPC1 = outdir / "PSPC1_RPKM.png"
    module.SUMMARY = outdir / "RPKM_profiles_summary.tsv"
    module.REGION_SPECS = (
        ("111", "111_MCM3_NONO_PSPC1", "MCM3-PSPC1-NONO co-binding", counts["111_MCM3_NONO_PSPC1"]),
        ("110", "110_MCM3_NONO", "MCM3-NONO co-binding", counts["110_MCM3_NONO"]),
        ("101", "101_MCM3_PSPC1", "MCM3-PSPC1 co-binding", counts["101_MCM3_PSPC1"]),
        ("100", "100_only_MCM3", "MCM3-specific binding", counts["100_only_MCM3"]),
    )
    module.FACTOR_SPECS = (
        ("MCM3", 0, module.OUTPUT),
        ("NONO", 1, module.OUTPUT_NONO),
        ("PSPC1", 2, module.OUTPUT_PSPC1),
    )
    def fitted_r_plot(values_tsv: Path, output_png: Path) -> None:
        values = pd.read_csv(values_tsv, sep="\t").groupby("group", observed=True)["value"]
        whiskers = []
        for _, group_values in values:
            finite = group_values[np.isfinite(group_values)].to_numpy(dtype=float)
            q1, q3 = np.quantile(finite, [0.25, 0.75])
            whiskers.append(float(finite[finite <= q3 + 1.5 * (q3 - q1)].max()))
        target = max(1.0, max(whiskers) * 1.08)
        magnitude = 10.0 ** math.floor(math.log10(target))
        module.Y_MAX = next(value * magnitude for value in (1, 1.25, 1.5, 2, 2.5, 3, 4, 5, 6, 8, 10) if value * magnitude >= target)
        labels = [
            f"MCM3-PSPC1-NONO co-binding\\n(n={counts['111_MCM3_NONO_PSPC1']:,})",
            f"MCM3-NONO co-binding\\n(n={counts['110_MCM3_NONO']:,})",
            f"MCM3-PSPC1 co-binding\\n(n={counts['101_MCM3_PSPC1']:,})",
            f"MCM3-specific binding\\n(n={counts['100_only_MCM3']:,})",
        ]
        levels = ["MCM3-PSPC1-NONO co-binding", "MCM3-NONO co-binding", "MCM3-PSPC1 co-binding", "MCM3-specific binding"]
        r_script = r'''
args <- commandArgs(trailingOnly = TRUE)
df <- read.delim(args[1], sep = "\t", header = TRUE, stringsAsFactors = FALSE)
levels <- strsplit(args[4], "\\|", fixed = FALSE)[[1]]
labels <- strsplit(args[5], "\\|", fixed = FALSE)[[1]]
df$group <- factor(df$group, levels = levels)
cols <- c("#cf111f", "#f2a010", "#4f7db1", "#777777")
png(args[2], width = 2400, height = 1816, res = 300, bg = "white")
par(mar = c(11, 7, 2, 2) + 0.1, cex.axis = 1.25, cex.lab = 1.45, font.axis = 2, font.lab = 2)
boxplot(value ~ group, data = df, col = cols, border = "#222222", lwd = 2.0,
        ylab = args[3], xlab = "", xaxt = "n", yaxt = "n", outline = FALSE,
        ylim = c(0, as.numeric(args[6])), whisklty = 1, staplelty = 1)
axis(1, at = seq_along(levels), labels = FALSE, tick = FALSE)
text(x = seq_along(levels), y = rep(par("usr")[3] - 0.39, length(levels)), labels = labels,
     srt = 35, xpd = TRUE, adj = 1, cex = 0.82, font = 2)
axis(2, at = seq(0, as.numeric(args[6]), by = 1), labels = seq(0, as.numeric(args[6]), by = 1), las = 1, cex.axis = 1.25, font = 2)
dev.off()
'''
        process = subprocess.run(
            [RSCRIPT, "-", str(values_tsv), str(output_png), module.SIGNAL_LABEL, "|".join(levels), "|".join(labels), str(module.Y_MAX)],
            input=r_script,
            text=True,
            capture_output=True,
            check=False,
        )
        if process.returncode:
            raise RuntimeError(process.stderr or process.stdout)

    module.r_plot = fitted_r_plot
    module.main()
    for name in ("MCM3_RPKM.png", "NONO_RPKM.png", "PSPC1_RPKM.png"):
        (outdir / name).replace(VISUALS / name)
    if module.SUMMARY.is_file():
        module.SUMMARY.replace(DATA / "rpkm_profiles" / module.SUMMARY.name)


def main() -> None:
    required = [GENE_BED, *TRACKS.values(), IGG]
    missing = [str(path) for path in required if not path.is_file() or path.stat().st_size == 0]
    if missing:
        raise FileNotFoundError("Missing rendering inputs:\n" + "\n".join(missing))
    for directory in (VISUALS, DATA / "promoter_profiles", DATA / "rpkm_profiles"):
        directory.mkdir(parents=True, exist_ok=True)
    sets = promoter_sets()
    render_promoter_venn(sets)
    render_promoter_profiles(sets)
    render_rpkm_profiles(sets)
    print(f"[DONE] CUT&RUN promoter and RPKM panels: {VISUALS}")


if __name__ == "__main__":
    main()
