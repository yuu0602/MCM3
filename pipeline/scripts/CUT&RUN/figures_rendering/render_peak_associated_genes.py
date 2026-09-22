#!/usr/bin/env python3
"""Render whole-genome nearest-TSS peak-associated gene figures."""

from bisect import bisect_left
import importlib.util
from pathlib import Path

import pandas as pd

PROJECT = Path(__file__).resolve().parent
RUN = PROJECT.parents[2]
CUTRUN = RUN / "cutrun_work"
DATA = CUTRUN / "data" / "figure_inputs"
VISUALS = CUTRUN / "visuals"
GTF = RUN / "reference" / "gencode.vM25.annotation.gtf"
FACTORS = ("MCM3", "NONO", "PSPC1")


def gene_tss() -> dict[str, list[tuple[int, str]]]:
    rows = []
    with GTF.open() as handle:
        for line in handle:
            fields = line.rstrip().split("\t")
            if len(fields) != 9 or fields[2] != "gene":
                continue
            attributes = fields[8]
            if 'gene_name "' not in attributes:
                continue
            symbol = attributes.split('gene_name "')[1].split('"')[0]
            start, end = int(fields[3]) - 1, int(fields[4])
            rows.append((fields[0], start if fields[6] == "+" else end, symbol))
    table = pd.DataFrame(rows, columns=["chrom", "tss", "gene"]).drop_duplicates()
    return {chrom: list(group[["tss", "gene"]].itertuples(index=False, name=None)) for chrom, group in table.sort_values(["chrom", "tss"]).groupby("chrom", sort=False)}


def nearest_gene(chrom: str, start: int, end: int, index: dict[str, list[tuple[int, str]]]) -> tuple[str, int] | None:
    entries = index.get(chrom)
    if not entries:
        return None
    center = (start + end) // 2
    positions = [entry[0] for entry in entries]
    point = bisect_left(positions, center)
    candidates = entries[max(0, point - 1):point + 1]
    tss, gene = min(candidates, key=lambda entry: (abs(entry[0] - center), entry[1]))
    return gene, abs(tss - center)


def peak_sets() -> dict[str, pd.DataFrame]:
    loci = pd.read_csv(CUTRUN / "data" / "PeakLoci" / "Venn_Peaks_loci.tsv", sep="\t")
    sets = {}
    for factor, bit in zip(FACTORS, (1, 2, 4)):
        sets[factor] = loci.loc[loci.region_mask.astype(int).map(lambda mask: bool(mask & bit)), ["chrom", "start", "end"]].copy()
    igg = pd.read_csv(CUTRUN / "data" / "PeakSets" / "Peaks_IgG.bed", sep="\t", header=None, names=["chrom", "start", "end"])
    sets["IgG"] = igg
    return sets


def load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot load {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def gene_annotation_renderer(output: Path):
    distribution = load_module(PROJECT / "render_peak_distribution.py", "peak_distribution")
    renderer = load_module(PROJECT / "peak_annotation.py", "peak_annotation")
    renderer.WORK_ROOT = CUTRUN
    renderer.OUTDIR = output / "annotation"
    renderer.TMPDIR = renderer.OUTDIR / "intermediate"
    renderer.PAPER_TMPDIR = renderer.OUTDIR / "paperfig_intermediate"
    renderer.GTF_CANDIDATES = [GTF]
    renderer.ASSAYS = list(FACTORS)
    renderer.PEAK_UNIT_LABEL = "peak-associated genes"
    renderer.USE_NT_KD_UNION = False
    renderer.BEDTOOLS = Path("/opt/anaconda3/envs/cutrun_env/bin/bedtools")
    renderer.ensure_dirs()
    return renderer, distribution.category_beds(renderer, CUTRUN / "data" / "PeakLoci" / "Venn_Peaks_MCM3.bed")


def annotate_unique_genes(renderer, categories: dict[str, Path], table: pd.DataFrame, factor: str) -> dict[str, int]:
    bed = renderer.TMPDIR / f"{factor}_unique_peak_associated_genes.bed"
    table[["chrom", "start", "end"]].to_csv(bed, sep="\t", header=False, index=False)
    midpoint = renderer.TMPDIR / f"{factor}_unique_peak_associated_genes_mid.bed"
    renderer.bed_to_peak_midpoints(bed, midpoint)
    return renderer.assign_categories(midpoint, categories)


def main() -> None:
    spec = importlib.util.spec_from_file_location("peak_venn", PROJECT / "peak_venn.py")
    if spec is None or spec.loader is None:
        raise RuntimeError("Cannot load peak_venn.py")
    venn = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(venn)

    output = DATA / "peak_associated_genes"
    output.mkdir(parents=True, exist_ok=True)
    VISUALS.mkdir(parents=True, exist_ok=True)
    index, assigned = gene_tss(), {}
    for factor, peaks in peak_sets().items():
        rows = []
        for row in peaks.itertuples(index=False):
            hit = nearest_gene(str(row.chrom), int(row.start), int(row.end), index)
            if hit:
                gene, distance = hit
                rows.append((str(row.chrom), int(row.start), int(row.end), gene, distance))
        table = pd.DataFrame(rows, columns=["chrom", "start", "end", "gene", "nearest_tss_distance_bp"])
        table = table.sort_values(["gene", "nearest_tss_distance_bp", "chrom", "start"]).drop_duplicates("gene")
        assigned[factor] = table
        table.to_csv(output / f"PeakAssociatedGenes_{factor}.tsv", sep="\t", index=False)
    renderer, categories = gene_annotation_renderer(output)
    counts = {
        factor: annotate_unique_genes(renderer, categories, table, factor)
        for factor, table in assigned.items()
    }
    for factor in assigned:
        renderer.save_single(
            VISUALS / f"Pie_PeakAssociatedGenes_{factor}.png",
            factor,
            counts[factor],
            len(assigned[factor]),
        )
    renderer.save_three_panel(
        VISUALS / "Pie_PeakAssociatedGenes.png",
        {factor: counts[factor] for factor in FACTORS},
        {factor: len(assigned[factor]) for factor in FACTORS},
    )
    annotation_rows = [
        {"factor": factor, "category": category, "n_unique_peak_associated_genes": value}
        for factor, result in counts.items()
        for category, value in result.items()
    ]
    pd.DataFrame(annotation_rows).to_csv(output / "PeakAssociatedGenes_annotation_counts.tsv", sep="\t", index=False)
    factor_sets = tuple(set(assigned[factor].gene) for factor in FACTORS)
    venn.plot_triple_venn(factor_sets, FACTORS, VISUALS / "Venn_PeakAssociatedGenes.png")
    venn.plot_triple_venn(factor_sets, FACTORS, VISUALS / "Venn_PeakAssociatedGenes_noNumbers.png", show_numbers=False, show_totals=False)
    print("[DONE] Whole-genome nearest-TSS peak-associated gene figures:", VISUALS)


if __name__ == "__main__":
    main()
