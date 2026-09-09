#!/usr/bin/env python3
"""Render CUT&RUN peak-overlap and gene-body profile figures."""

from __future__ import annotations

import importlib.util
from pathlib import Path
import subprocess

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


PROJECT = Path(__file__).resolve().parent
RUN = PROJECT.parents[2]
CUTRUN = RUN / "cutrun_work"
DATA = CUTRUN / "data" / "figure_inputs"
VISUALS = CUTRUN / "visuals"
TAG = "p1e4_q5e2_fe3_union4of6"
GTF = RUN / "reference" / "gencode.vM25.annotation.gtf"
GENE_BED = CUTRUN / "data" / "GeneBodies_M25.bed6"
FACTORS = ("MCM3", "NONO", "PSPC1")
COLORS = {"MCM3": "#C4161C", "NONO": "#F39C12", "PSPC1": "#4C78A8", "IgG": "#737B8C"}
TRACKS = {
    factor: CUTRUN / "03_bigwig" / "IGV_representation" / f"{factor}_six_library_batchBalanced_mean.bw"
    for factor in FACTORS
}
PEAK_LOCI = CUTRUN / "data" / "PeakLoci" / "Venn_Peaks_loci.tsv"
PROMOTERS = CUTRUN / "data" / "Promoters_M25_TSSplusminus1kb.bed"
PROMOTER_OVERLAP_BP = 250
OUT = DATA / "additional_panels"
COMPUTE_MATRIX = Path("/opt/anaconda3/envs/cutrun_env/bin/computeMatrix")


def load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot import {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def save(fig: plt.Figure, name: str) -> None:
    VISUALS.mkdir(parents=True, exist_ok=True)
    fig.savefig(VISUALS / name, dpi=300, facecolor="white")
    plt.close(fig)


def promoter_sets() -> dict[str, set[str]]:
    root = DATA / "promoter_gene_venn"
    return {
        factor: set(pd.read_csv(root / f"Venn_PromoterGenes_{factor}.tsv", sep="\t", dtype=str)["gene"].dropna())
        for factor in FACTORS
    }


def draw_reference_venn() -> None:
    counts = pd.read_csv(DATA / "Venn_Peaks_counts.tsv", sep="\t")
    if "peak_loci" not in counts:
        raise RuntimeError("Peak Venn input must contain nonredundant peak_loci counts")
    by_mask = dict(zip(counts["region_mask"].astype(int), counts["peak_loci"].astype(int)))
    def tokens(prefix: str, size: int) -> set[str]:
        return {f"{prefix}{index}" for index in range(size)}
    only_mcm3 = tokens("m", by_mask[1]); only_nono = tokens("n", by_mask[2]); only_pspc1 = tokens("p", by_mask[4])
    mn = tokens("mn", by_mask[3]); mp = tokens("mp", by_mask[5]); npair = tokens("np", by_mask[6]); triple = tokens("all", by_mask[7])
    sets = {
        "MCM3": only_mcm3 | mn | mp | triple,
        "NONO": only_nono | mn | npair | triple,
        "PSPC1": only_pspc1 | mp | npair | triple,
    }
    venn = load_module(PROJECT / "peak_venn.py", "peak_venn")
    values = tuple(sets[factor] for factor in FACTORS)
    venn.plot_triple_venn(
        values, FACTORS, VISUALS / "FigS2e_peak_overlap_whole_genome.png",
    )
    venn.plot_triple_venn(
        values, FACTORS, VISUALS / "FigS2e_peak_overlap_whole_genome_no_numbers.png",
        show_numbers=False, show_totals=False,
    )



def matrix(path: Path, genes_bed: Path, tracks: list[Path], mode: str) -> tuple[list[str], np.ndarray]:
    if not COMPUTE_MATRIX.is_file():
        raise FileNotFoundError(COMPUTE_MATRIX)
    command = [str(COMPUTE_MATRIX), mode, "-S", *(str(track) for track in tracks), "-R", str(genes_bed)]
    if mode == "scale-regions":
        command.extend(("--beforeRegionStartLength", "3000", "--regionBodyLength", "2000", "--afterRegionStartLength", "1000"))
    else:
        command.extend(("--referencePoint", "TSS", "-b", "3000", "-a", "3000"))
    command.extend(("--binSize", "25", "--missingDataAsZero", "-p", "2", "-o", str(path)))
    subprocess.run(command, check=True)
    parser = load_module(PROJECT / "annotation.py", "annotation")
    return parser.read_matrix(path, len(tracks))



def peak_locus_groups() -> dict[int, pd.DataFrame]:
    """Load disjoint peak-locus Venn regions for the peak-centered panels."""
    table = pd.read_csv(PEAK_LOCI, sep="\t")
    required = {"region_mask", "chrom", "start", "end"}
    missing = required - set(table.columns)
    if missing:
        raise RuntimeError(f"Peak locus table lacks columns: {', '.join(sorted(missing))}")
    table["region_mask"] = pd.to_numeric(table["region_mask"], errors="raise").astype(int)
    table["start"] = pd.to_numeric(table["start"], errors="raise").astype(int)
    table["end"] = pd.to_numeric(table["end"], errors="raise").astype(int)
    columns = [column for column in ("locus_id", "chrom", "start", "end") if column in table.columns]
    return {mask: table.loc[table.region_mask.eq(mask), columns].copy() for mask in (1, 3, 5, 7)}



def peak_promoter_associations(groups: dict[int, pd.DataFrame]) -> pd.DataFrame:
    """Map Venn peak loci to M25 promoters using the final >=250-bp rule."""
    promoters = pd.read_csv(
        PROMOTERS, sep="\t", header=None, usecols=[0, 1, 2, 3],
        names=["chrom", "start", "end", "gene"], dtype={"chrom": str, "gene": str},
    )
    promoters["start"] = pd.to_numeric(promoters["start"], errors="raise").astype(int)
    promoters["end"] = pd.to_numeric(promoters["end"], errors="raise").astype(int)
    by_chrom = {
        chrom: table.sort_values(["start", "end", "gene"], kind="stable").reset_index(drop=True)
        for chrom, table in promoters.groupby("chrom", sort=False)
    }
    rows: list[dict[str, object]] = []
    for mask, loci in groups.items():
        for chrom, peak_table in loci.groupby("chrom", sort=False):
            promoter_table = by_chrom.get(str(chrom))
            if promoter_table is None:
                continue
            active: list[tuple[int, int, str]] = []
            next_promoter = 0
            for locus in peak_table.sort_values(["start", "end"], kind="stable").itertuples(index=False):
                peak_start, peak_end = int(locus.start), int(locus.end)
                while next_promoter < len(promoter_table) and int(promoter_table.iloc[next_promoter].start) <= peak_end - PROMOTER_OVERLAP_BP:
                    candidate = promoter_table.iloc[next_promoter]
                    active.append((int(candidate.start), int(candidate.end), str(candidate.gene)))
                    next_promoter += 1
                active = [candidate for candidate in active if candidate[1] >= peak_start + PROMOTER_OVERLAP_BP]
                for promoter_start, promoter_end, gene in active:
                    overlap = min(peak_end, promoter_end) - max(peak_start, promoter_start)
                    if overlap >= PROMOTER_OVERLAP_BP:
                        rows.append({
                            "region_mask": mask,
                            "locus_id": getattr(locus, "locus_id", f"{chrom}:{peak_start}-{peak_end}"),
                            "chrom": str(chrom),
                            "peak_start": peak_start,
                            "peak_end": peak_end,
                            "gene": gene,
                            "promoter_overlap_bp": overlap,
                        })
    return pd.DataFrame(rows, columns=["region_mask", "locus_id", "chrom", "peak_start", "peak_end", "gene", "promoter_overlap_bp"])


def write_peak_profile_summary(rows: list[dict[str, object]]) -> None:
    """Append peak-associated gene profile provenance to the centered-locus table."""
    current = CUTRUN / "data" / "Peak_profiles_sets.tsv"
    previous = pd.read_csv(current, sep="\t") if current.is_file() else pd.DataFrame()
    replacement = pd.DataFrame(rows)
    if not previous.empty and "panel" in previous:
        previous = previous.loc[~previous.panel.isin(set(replacement.panel))]
    summary = pd.concat((previous, replacement), ignore_index=True, sort=False)
    summary.to_csv(current, sep="\t", index=False)
    summary.to_csv(OUT / "Peak_profiles_sets.tsv", sep="\t", index=False)


def write_profile_gene_bed(genes: set[str], path: Path) -> int:
    """Write one canonical gene-body record per promoter-associated gene."""
    bodies = pd.read_csv(
        GENE_BED, sep="\t", header=None,
        names=["chrom", "start", "end", "gene", "score", "strand"], dtype={"chrom": str, "gene": str},
    )
    selected = bodies.loc[bodies.gene.isin(genes)].drop_duplicates("gene").sort_values(["chrom", "start", "end", "gene"], kind="stable")
    missing = genes - set(selected.gene)
    if missing:
        raise RuntimeError(f"{len(missing)} promoter-associated genes are missing from the M25 gene-body table")
    selected.to_csv(path, sep="\t", header=False, index=False)
    return len(selected)


def canonical_peak_promoter_associations(associations: pd.DataFrame) -> pd.DataFrame:
    """Assign each promoter-overlapping peak locus one deterministic promoter."""
    return (
        associations.sort_values(
            ["region_mask", "locus_id", "promoter_overlap_bp", "gene"],
            ascending=[True, True, False, True], kind="stable",
        )
        .drop_duplicates(["region_mask", "locus_id"], keep="first")
        .reset_index(drop=True)
    )


def write_peak_association_gene_bed(associations: pd.DataFrame, path: Path) -> int:
    """Write one strand-aware TSS/TES region per canonical peak-promoter locus."""
    bodies = pd.read_csv(
        GENE_BED, sep="\t", header=None,
        names=["chrom", "start", "end", "gene", "score", "strand"], dtype={"chrom": str, "gene": str},
    ).drop_duplicates("gene")
    selected = associations.loc[:, ["locus_id", "gene"]].merge(bodies, on="gene", how="left", validate="many_to_one", sort=False)
    if selected[["chrom", "start", "end", "strand"]].isna().any(axis=None):
        missing = selected.loc[selected["chrom"].isna(), "gene"].drop_duplicates().tolist()
        raise RuntimeError(f"{len(missing)} peak-associated promoters are missing M25 gene-body records")
    selected.loc[:, ["chrom", "start", "end", "locus_id", "score", "strand"]].to_csv(path, sep="\t", header=False, index=False)
    return len(selected)


def render_peak_associated_gene_profiles() -> None:
    """Render TSS/TES profiles for promoter genes associated with peak loci."""
    OUT.mkdir(parents=True, exist_ok=True)
    groups = peak_locus_groups()
    associations = peak_promoter_associations(groups)
    associations.to_csv(CUTRUN / "data" / "Peak_profiles_promoter_overlaps.tsv", sep="\t", index=False)
    associations.to_csv(OUT / "Peak_profiles_promoter_overlaps.tsv", sep="\t", index=False)
    canonical = canonical_peak_promoter_associations(associations)
    canonical.to_csv(CUTRUN / "data" / "Peak_profiles_promoter_assignments.tsv", sep="\t", index=False)
    canonical.to_csv(OUT / "Peak_profiles_promoter_assignments.tsv", sep="\t", index=False)
    if any(canonical.loc[canonical.region_mask.eq(mask)].empty for mask in groups):
        empty = [str(mask) for mask in groups if canonical.loc[canonical.region_mask.eq(mask)].empty]
        raise RuntimeError("No promoter-associated genes were recovered for peak Venn region(s): " + ", ".join(empty))

    # Match computeMatrix exactly: 3 kb upstream / 25 bp = 120 bins,
    # 2 kb scaled body / 25 bp = 80 bins, and 1 kb downstream = 40 bins.
    x_body = np.concatenate((
        np.linspace(-3.0, 0.0, 120, endpoint=False),
        np.linspace(0.0, 2.0, 80, endpoint=False),
        np.linspace(2.0, 3.0, 40, endpoint=False),
    ))
    summary_rows: list[dict[str, object]] = []

    group_specs = (
        ("Triple co-localized", 7, COLORS["MCM3"]),
        ("MCM3-NONO", 3, COLORS["NONO"]),
        ("MCM3-PSPC1", 5, COLORS["PSPC1"]),
        ("MCM3-specific", 1, "#777777"),
    )
    fig = plt.figure(figsize=(7.4, 5.6), dpi=300, facecolor="white")
    axis = fig.add_axes([.11, .20, .57, .62])
    profiles: dict[str, np.ndarray] = {}
    ymax = 0.0
    for label, mask, color in group_specs:
        peak_subset = canonical.loc[canonical.region_mask.eq(mask)].copy()
        peak_bed = OUT / f"MCM3_Peak_profiles_region_{mask:03b}_regions.bed"
        write_peak_association_gene_bed(peak_subset, peak_bed)
        _, values = matrix(OUT / f"MCM3_Peak_profiles_region_{mask:03b}_matrix.gz", peak_bed, [TRACKS["MCM3"]], "scale-regions")
        profile = values[:, 0, :].mean(axis=0)
        retained = values.shape[0]
        if retained != len(peak_subset):
            raise RuntimeError(f"{label}: retained {retained} loci but expected {len(peak_subset)} canonical promoter-overlapping peaks")
        profiles[label] = profile
        ymax = max(ymax, float(np.max(profile)))
        axis.plot(x_body, profile, lw=2.6, color=color, label=f"{label}\npeaks")
        summary_rows.append({
            "panel": "MCM3_Peak_profiles", "group": label, "region_mask": mask,
            "unit": "promoter-overlapping peak loci", "n_input_peak_loci": len(groups[mask]),
            "n_profile_loci": retained, "promoter_overlap_rule": f">={PROMOTER_OVERLAP_BP} bp",
            "reference_point": "TSS/TES scaled gene body",
        })
    axis.axvline(0.0, color="#222222", lw=1.1, ls="--"); axis.axvline(2.0, color="#222222", lw=1.1, ls="--")
    axis.set_xlim(-3.0, 3.0); axis.set_ylim(0.0, ymax * 1.05)
    axis.set_xticks((-3.0, 0.0, 2.0, 3.0), ("-3.0", "TSS", "TES", "3.0"))
    axis.set_xlabel("Relative distance (kb)", fontsize=16, fontweight="bold"); axis.set_ylabel("Coverage", fontsize=16, fontweight="bold")
    axis.tick_params(labelsize=14)
    for tick_label in (*axis.get_xticklabels(), *axis.get_yticklabels()): tick_label.set_fontweight("bold")
    legend = axis.legend(frameon=False, loc="center left", bbox_to_anchor=(1.03, .5), fontsize=8.5)
    for label in legend.get_texts(): label.set_fontweight("bold")
    save(fig, "MCM3_Peak_profiles.png")
    pd.DataFrame({"relative_position_kb": x_body, **profiles}).to_csv(
        OUT / "MCM3_Peak_profiles_data.tsv", sep="\t", index=False
    )

    triple_peaks = canonical.loc[canonical.region_mask.eq(7)].copy()
    triple_bed = OUT / "MCM3_NONO_PSPC1_peak_profiles_regions.bed"
    write_peak_association_gene_bed(triple_peaks, triple_bed)
    _, matrix_values = matrix(OUT / "MCM3_NONO_PSPC1_peak_profiles_matrix.gz", triple_bed, [TRACKS[factor] for factor in FACTORS], "scale-regions")
    values = {factor: matrix_values[:, index, :].mean(axis=0) for index, factor in enumerate(FACTORS)}
    retained = matrix_values.shape[0]
    if retained != len(triple_peaks):
        raise RuntimeError(f"Triple co-localized: retained {retained} loci but expected {len(triple_peaks)} canonical promoter-overlapping peaks")
    fig = plt.figure(figsize=(7.4, 5.6), dpi=300, facecolor="white")
    axis = fig.add_axes([.11, .20, .64, .62])
    for factor in FACTORS:
        axis.plot(x_body, values[factor], lw=2.7, color=COLORS[factor], label=factor)
    axis.axvline(0.0, color="#222222", lw=1.1, ls="--"); axis.axvline(2.0, color="#222222", lw=1.1, ls="--")
    axis.set_xlim(-3.0, 3.0); axis.set_ylim(0.0, max(float(np.max(profile)) for profile in values.values()) * 1.05)
    axis.set_xticks((-3.0, 0.0, 2.0, 3.0), ("-3.0", "TSS", "TES", "3.0"))
    axis.set_xlabel("Relative distance (kb)", fontsize=16, fontweight="bold"); axis.set_ylabel("Coverage", fontsize=16, fontweight="bold")
    axis.tick_params(labelsize=14)
    for tick_label in (*axis.get_xticklabels(), *axis.get_yticklabels()): tick_label.set_fontweight("bold")
    legend = axis.legend(frameon=False, loc="center left", bbox_to_anchor=(1.04, .5), fontsize=12)
    for label in legend.get_texts(): label.set_fontweight("bold")
    save(fig, "MCM3_NONO_PSPC1_peak_profiles.png")
    pd.DataFrame({"relative_position_kb": x_body, **values}).to_csv(
        OUT / "MCM3_NONO_PSPC1_peak_profiles_data.tsv", sep="\t", index=False
    )
    summary_rows.append({
        "panel": "MCM3_NONO_PSPC1_peak_profiles", "group": "MCM3+NONO+PSPC1", "region_mask": 7,
        "unit": "promoter-overlapping peak loci", "n_input_peak_loci": len(groups[7]),
        "n_profile_loci": retained, "promoter_overlap_rule": f">={PROMOTER_OVERLAP_BP} bp",
        "reference_point": "TSS/TES scaled gene body",
    })
    write_peak_profile_summary(summary_rows)



def style_axis(axis: plt.Axes, tick_size: int) -> None:
    axis.tick_params(axis="both", width=1.6, length=5, labelsize=tick_size)
    for label in (*axis.get_xticklabels(), *axis.get_yticklabels()):
        label.set_fontweight("bold")
    for spine in axis.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(1.6)


def render_promoter_gene_profile() -> None:
    sets = promoter_sets()
    genes = set().union(*sets.values())
    region_bed = OUT / "Metaprofile_PromoterGenes_regions.bed"
    annotation = load_module(PROJECT / "annotation.py", "promoter_annotation")
    n_genes = annotation.write_gene_bed(GTF, genes, region_bed)
    _, values = matrix(
        OUT / "Metaprofile_PromoterGenes_matrix.gz",
        region_bed,
        [TRACKS[factor] for factor in FACTORS],
        "scale-regions",
    )
    profiles = {factor: values[:, index, :].mean(axis=0) for index, factor in enumerate(FACTORS)}
    x = np.concatenate((
        np.linspace(-3.0, 0.0, 120, endpoint=False),
        np.linspace(0.0, 2.0, 80, endpoint=False),
        np.linspace(2.0, 3.0, 40, endpoint=False),
    ))
    pd.DataFrame({"relative_position_kb": x, **profiles}).to_csv(
        OUT / "Metaprofile_PromoterGenes_data.tsv", sep="\t", index=False
    )
    fig = plt.figure(figsize=(9.5, 9.5 * 5.6 / 7.4), facecolor="white")
    axis = fig.add_subplot(111)
    for factor in FACTORS:
        axis.plot(x, profiles[factor], lw=2.7, label=factor, color=COLORS[factor])
    axis.axvline(0.0, ls="--", lw=1.6, color="black")
    axis.axvline(2.0, ls="--", lw=1.6, color="black")
    axis.set_xlim(-3.0, 3.0)
    axis.set_ylim(0.0, max(float(np.max(profile)) for profile in profiles.values()) * 1.05)
    axis.set_xticks((-3.0, 0.0, 2.0, 3.0), ("-3.0", "TSS", "TES", "3.0"))
    axis.set_title("CUT&RUN promoter-bound genes", fontsize=14, fontweight="bold", pad=18)
    axis.text(
        .5, 1.012,
        f"p ≤ 1e-4; total support ≥4/6; overlap ≥250 bp; union n={n_genes:,}",
        transform=axis.transAxes, ha="center", va="bottom", fontsize=9, fontweight="bold",
    )
    axis.set_xlabel("Relative position (kb)", fontsize=15, fontweight="bold")
    axis.set_ylabel("Coverage", fontsize=15, fontweight="bold")
    style_axis(axis, 13)
    legend = axis.legend(frameon=False, fontsize=12, loc="center left", bbox_to_anchor=(1.02, .5))
    for label in legend.get_texts():
        label.set_fontweight("bold")
    fig.subplots_adjust(left=.13, right=.76, bottom=.16, top=.79)
    save(fig, "FigX_TSS_to_TES_metaprofile__NT_KD_UNION__minus3_to_plus3_TES2p8.png")


def main() -> None:
    required = [
        GTF, GENE_BED, COMPUTE_MATRIX, *TRACKS.values(),
        DATA / "Venn_Peaks_counts.tsv", PEAK_LOCI, PROMOTERS,
    ]
    missing = [str(path) for path in required if not path.is_file() or path.stat().st_size == 0]
    if missing:
        raise FileNotFoundError("Missing reference-layout input(s):\n" + "\n".join(missing))
    draw_reference_venn()
    render_promoter_gene_profile()
    render_peak_associated_gene_profiles()
    print(f"[DONE] CUT&RUN overlap and metaprofile figures: {VISUALS}")


if __name__ == "__main__":
    main()
