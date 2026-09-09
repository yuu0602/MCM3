#!/usr/bin/env python3
"""Render gene-body profiles for promoter-gene Venn regions."""

from __future__ import annotations

import gzip
import argparse
import subprocess
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.ticker import MaxNLocator


WORK_ROOT = Path("cutrun_work_20191012")
FIGURE_DATA = WORK_ROOT / "data" / "figure_inputs"
OUTDIR = WORK_ROOT / "visuals"
TMPDIR = FIGURE_DATA / "promoter_profiles" / "intermediate"
REGION_TABLE = FIGURE_DATA / "promoter_gene_venn" / "Venn_PromoterGenes_regions.tsv"
GENE_BED = WORK_ROOT / "data" / "GeneBodies_M25.bed6"
COMPUTE_MATRIX = Path("/opt/anaconda3/envs/cutrun_env/bin/computeMatrix")
LOCAL_BIGWIG_CACHE = Path("/tmp/mcm3_cutrun_bigwig_cache_20191012")
LOCAL_GENE_BED = LOCAL_BIGWIG_CACHE / "genes.norm.sorted.bed6"
MATRIX_TMPDIR = Path("/tmp/mcm3_cutrun_promoter_profile_matrices")

ASSAYS = ("MCM3", "NONO", "PSPC1")
BIGWIGS = {
    "MCM3": WORK_ROOT / "03_bigwig/NT_rep1_MCM3/NT_rep1_MCM3.bw",
    "NONO": WORK_ROOT / "03_bigwig/NT_rep1_NONO/NT_rep1_NONO.bw",
    "PSPC1": WORK_ROOT / "03_bigwig/NT_rep1_PSPC1/NT_rep1_PSPC1.bw",
}
REGIONS = {
    "111": {
        "data_stem": "Profile_MCM3_NONO_PSPC1",
        "table_region": "111_MCM3_NONO_PSPC1",
        "label": "MCM3-NONO-PSPC1 co-bound promoter genes",
        "display_title": "MCM3-NONO-PSPC1 Co-localized\ngenes",
        "tracks": ("MCM3", "NONO", "PSPC1"),
        "expected_n": 3473,
    },
    "110": {
        "data_stem": "Profile_MCM3_NONO",
        "table_region": "110_MCM3_NONO",
        "label": "MCM3-NONO-only promoter genes",
        "display_title": "MCM3-NONO-only\npromoter genes",
        "tracks": ("MCM3", "NONO"),
        "expected_n": 177,
    },
    "101": {
        "data_stem": "Profile_MCM3_PSPC1",
        "table_region": "101_MCM3_PSPC1",
        "label": "MCM3-PSPC1-only promoter genes",
        "display_title": "MCM3-PSPC1-only\npromoter genes",
        "tracks": ("MCM3", "PSPC1"),
        "expected_n": 2307,
    },
}

UPSTREAM_BP = 3000
XMAX_BP = 3000
BODY_BP = 2000
DOWNSTREAM_BP = XMAX_BP - BODY_BP
BIN_SIZE = 25
DPI = 300
UPPER_TAIL_TRIM = 0.01
DISPLAY_SCALE = 1000.0 / BIN_SIZE
SIGNAL_LABEL = "Coverage (RPKM)"
# Match the established promoter-profile reference figure.
COLORS = {"MCM3": "#A80F14", "NONO": "#F39C12", "PSPC1": "#1F4AA8"}


def fail(message: str) -> None:
    raise RuntimeError(message)


def run(command: list[str]) -> None:
    print("[CMD] " + " ".join(command))
    subprocess.run(command, check=True)


def load_region_genes() -> dict[str, set[str]]:
    table = pd.read_csv(REGION_TABLE, sep="\t", dtype=str)
    if set(table.columns) != {"region", "gene"}:
        fail(f"Unexpected region-table columns: {list(table.columns)}")

    region_genes: dict[str, set[str]] = {}
    for code, config in REGIONS.items():
        genes = set(table.loc[table["region"] == config["table_region"], "gene"].dropna())
        if len(genes) != config["expected_n"]:
            fail(
                f"Region {code} contains {len(genes):,} genes; "
                f"expected {config['expected_n']:,} from the p1e3/min2of4 Venn."
            )
        region_genes[code] = genes

    return region_genes


def resolve_profile_bigwigs() -> dict[str, Path]:
    """Prefer byte-matched local copies when manually staged for external-drive I/O."""
    cached = {assay: LOCAL_BIGWIG_CACHE / source.name for assay, source in BIGWIGS.items()}
    if all(path.is_file() and path.stat().st_size == BIGWIGS[assay].stat().st_size for assay, path in cached.items()):
        print(f"[OK] Using byte-matched local bigWig cache: {LOCAL_BIGWIG_CACHE}")
        return cached
    print("[OK] Using source bigWigs from the CUT&RUN workspace.")
    return BIGWIGS


def resolve_gene_bed() -> Path:
    if LOCAL_GENE_BED.is_file() and LOCAL_GENE_BED.stat().st_size == GENE_BED.stat().st_size:
        print(f"[OK] Using byte-matched local gene-body BED: {LOCAL_GENE_BED}")
        return LOCAL_GENE_BED
    print(f"[OK] Using source gene-body BED: {GENE_BED}")
    return GENE_BED


def write_gene_beds(region_genes: dict[str, set[str]], annotation_bed: Path) -> dict[str, Path]:
    """Subset the normalized GTF-derived gene-body BED used by earlier profiles."""
    wanted = set().union(*region_genes.values())
    gene_bed = pd.read_csv(
        annotation_bed,
        sep="\t",
        header=None,
        names=["chrom", "start", "end", "gene", "score", "strand"],
        dtype={"chrom": str, "start": int, "end": int, "gene": str, "score": str, "strand": str},
    )
    found = set(gene_bed.loc[gene_bed["gene"].isin(wanted), "gene"])
    missing = wanted - found
    if missing:
        preview = ", ".join(sorted(missing)[:20])
        fail(f"{len(missing):,} Venn genes are absent from the GTF, including: {preview}")

    out_beds: dict[str, Path] = {}
    for code, genes in region_genes.items():
        expected_n = REGIONS[code]["expected_n"]
        region_bed = gene_bed.loc[gene_bed["gene"].isin(genes)].copy()
        n_unique = region_bed["gene"].nunique()
        if n_unique != expected_n:
            fail(f"Region {code} yielded {n_unique:,} unique gene records; expected {expected_n:,}.")
        region_bed = region_bed.sort_values(["chrom", "start", "end", "gene"], kind="stable")
        out_bed = TMPDIR / f"{REGIONS[code]['data_stem']}_regions.bed"
        region_bed.to_csv(out_bed, sep="\t", header=False, index=False)
        out_beds[code] = out_bed
    return out_beds


def compute_matrix(region_bed: Path, out_matrix: Path, profile_bigwigs: dict[str, Path]) -> None:
    run(
        [
            str(COMPUTE_MATRIX),
            "scale-regions",
            "-S",
            *[str(profile_bigwigs[assay]) for assay in ASSAYS],
            "-R",
            str(region_bed),
            "--beforeRegionStartLength",
            str(UPSTREAM_BP),
            "--regionBodyLength",
            str(BODY_BP),
            "--afterRegionStartLength",
            str(DOWNSTREAM_BP),
            "--binSize",
            str(BIN_SIZE),
            "--missingDataAsZero",
            # macOS spawned workers can block indefinitely while reading these
            # external-volume bigWigs; one worker is slower but completes reliably.
            "-p",
            "1",
            "-o",
            str(out_matrix),
        ]
    )


def read_trimmed_mean_profiles(matrix_gz: Path) -> dict[str, np.ndarray]:
    rows: list[np.ndarray] = []
    with gzip.open(matrix_gz, "rt") as handle:
        for line in handle:
            if line.startswith(("@", "#")) or not line.strip():
                continue
            values = np.array(
                [float(value) if value != "nan" else np.nan for value in line.rstrip("\n").split("\t")[6:]],
                dtype=float,
            )
            rows.append(values)

    if not rows:
        fail(f"computeMatrix produced no data rows: {matrix_gz}")
    matrix = np.vstack(rows)
    if matrix.shape[1] % len(ASSAYS) != 0:
        fail(f"Unexpected matrix width {matrix.shape[1]} for {len(ASSAYS)} tracks: {matrix_gz}")
    n_bins = matrix.shape[1] // len(ASSAYS)
    matrix = np.nan_to_num(matrix.reshape(matrix.shape[0], len(ASSAYS), n_bins), nan=0.0)
    # Trim only the upper 1% per bin so isolated high-coverage loci do not
    # dominate the mean while all Venn genes remain represented.
    # per bin so all Venn genes remain represented without those artifacts dominating.
    n_keep = int(np.floor(matrix.shape[0] * (1.0 - UPPER_TAIL_TRIM)))
    if n_keep < 1:
        fail(f"Too few regions for {UPPER_TAIL_TRIM:.0%} upper-tail trimming: {matrix.shape[0]}")
    trimmed_mean = np.sort(matrix, axis=0)[:n_keep, :, :].mean(axis=0) * DISPLAY_SCALE
    return {assay: trimmed_mean[index, :] for index, assay in enumerate(ASSAYS)}


def plot_profile(code: str, profiles: dict[str, np.ndarray]) -> Path:
    n_bins = len(next(iter(profiles.values())))
    n_up = int(round(UPSTREAM_BP / BIN_SIZE))
    n_body = int(round(BODY_BP / BIN_SIZE))
    n_down = int(round(DOWNSTREAM_BP / BIN_SIZE))
    if n_up + n_body + n_down != n_bins:
        fail(f"Unexpected bin count: {n_bins}; expected {n_up + n_body + n_down}.")
    # Reference display coordinates: TSS and TES are internal axis landmarks.
    x_kb = np.concatenate(
        [
            np.linspace(-3.0, -1.5, n_up, endpoint=False),
            np.linspace(-1.5, 1.5, n_body, endpoint=False),
            np.linspace(1.5, 3.0, n_down, endpoint=False),
        ]
    )
    fig = plt.figure(figsize=(7.4, 5.6), dpi=DPI)
    axis = fig.add_axes([0.11, 0.20, 0.64, 0.62])

    for assay in REGIONS[code]["tracks"]:
        axis.plot(x_kb, profiles[assay], lw=2.6, label=assay, color=COLORS[assay])
    axis.set_xlim(-UPSTREAM_BP / 1000.0, XMAX_BP / 1000.0)
    maximum = max(float(np.max(profiles[assay])) for assay in REGIONS[code]["tracks"])
    axis.set_ylim(0.0, maximum * 1.05)
    axis.axvline(-1.5, color="#222222", lw=1.1, ls="--")
    axis.axvline(1.5, color="#222222", lw=1.1, ls="--")
    axis.set_xticks([-3.0, -1.5, 1.5, 3.0])
    axis.set_xticklabels(["-3.0", "TSS", "TES", "3.0"])
    axis.set_xlabel("Relative distance (kb)", fontsize=16, fontweight="bold")
    axis.set_ylabel(SIGNAL_LABEL, fontsize=16, fontweight="bold")
    axis.yaxis.set_major_locator(MaxNLocator(integer=True, nbins=6))
    axis.tick_params(labelsize=14)
    legend = axis.legend(
        frameon=False,
        fontsize=15,
        loc="center left",
        bbox_to_anchor=(1.04, 0.5),
        borderaxespad=0.0,
    )
    for text in legend.get_texts():
        text.set_fontweight("bold")
    title_text = f"{REGIONS[code]['display_title']} (n={REGIONS[code]['expected_n']:,})"
    fig.suptitle(title_text, fontsize=22, fontweight="bold", x=0.5, y=0.97)
    for label in (*axis.get_xticklabels(), *axis.get_yticklabels()):
        label.set_fontweight("bold")

    out_png = OUTDIR / f"FigX_promoter_only_profile__region_{code}.png"
    fig.savefig(out_png, dpi=DPI, facecolor="white")
    plt.close(fig)
    return out_png


def main(reuse_matrices: bool = False) -> None:
    required = [REGION_TABLE, GENE_BED, COMPUTE_MATRIX, *BIGWIGS.values()]
    missing = [str(path) for path in required if not path.is_file() or path.stat().st_size == 0]
    if missing:
        fail("Missing required input(s):\n" + "\n".join(missing))

    OUTDIR.mkdir(parents=True, exist_ok=True)
    TMPDIR.mkdir(parents=True, exist_ok=True)
    MATRIX_TMPDIR.mkdir(parents=True, exist_ok=True)
    region_genes = load_region_genes()
    annotation_bed = resolve_gene_bed()
    gene_beds = write_gene_beds(region_genes, annotation_bed)
    profile_bigwigs = resolve_profile_bigwigs()

    summary_rows = []
    for code in REGIONS:
        matrix_path = MATRIX_TMPDIR / f"{REGIONS[code]['data_stem']}_matrix.gz"
        if reuse_matrices and matrix_path.is_file() and matrix_path.stat().st_size > 0:
            print(f"[OK] Reusing validated matrix: {matrix_path}")
        else:
            compute_matrix(gene_beds[code], matrix_path, profile_bigwigs)
        out_png = plot_profile(code, read_trimmed_mean_profiles(matrix_path))
        summary_rows.append(
            {
                "region_code": code,
                "venn_region": REGIONS[code]["table_region"],
                "n_genes": REGIONS[code]["expected_n"],
                "gene_bed": str(gene_beds[code]),
                "n_gene_body_records": sum(1 for _ in gene_beds[code].open()),
                "signal_aggregation": "upper_1pct_trimmed_mean",
                "display_scale": DISPLAY_SCALE,
                "signal_label": SIGNAL_LABEL,
                "profile_png": str(out_png),
            }
        )
        print(f"[OK] Wrote {out_png}")

    pd.DataFrame(summary_rows).to_csv(
        OUTDIR / "Profile_PromoterGenes_summary.tsv", sep="\t", index=False
    )
    print(f"[DONE] Wrote Venn-consistent promoter-gene profiles to {OUTDIR}")


if __name__ == "__main__":
    try:
        parser = argparse.ArgumentParser()
        parser.add_argument(
            "--reuse-matrices",
            action="store_true",
            help="Regenerate figures from current matrices only when their scale-region parameters are unchanged.",
        )
        main(reuse_matrices=parser.parse_args().reuse_matrices)
    except (RuntimeError, subprocess.CalledProcessError) as error:
        print(f"[ERROR] {error}", file=sys.stderr)
        raise SystemExit(1)
