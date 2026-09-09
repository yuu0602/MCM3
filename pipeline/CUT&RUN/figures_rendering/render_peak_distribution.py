#!/usr/bin/env python3
"""Annotate peak loci and render peak-distribution pies."""

from __future__ import annotations

import importlib.util
from pathlib import Path

import pandas as pd


PROJECT = Path(__file__).resolve().parent
RUN = PROJECT.parents[2]
CUTRUN = RUN / "cutrun_work"
DATA = CUTRUN / "data" / "figure_inputs"
VISUALS = CUTRUN / "visuals"
PEAKS = DATA / "batch_stratified_peaks"
STAGE = DATA / "peak_distribution" / "intermediate"
GTF = RUN / "reference" / "gencode.vM25.annotation.gtf"
FACTORS = ("MCM3", "NONO", "PSPC1")
IGG_CALLS = (
    CUTRUN / "04_peaks" / "2020" / "rep1" / "IgG" / "IgG-1_full_peaks.narrowPeak",
    CUTRUN / "04_peaks" / "2020" / "rep2" / "IgG" / "IgG-2_full_peaks.narrowPeak",
)
# A raw rebuild writes the strict two-IgG union directly. The accepted path
# leaves this unset and retains the original per-call union procedure.
IGG_BED: Path | None = None
# When supplied by the final workflow, these are the canonical merged loci
# used by the whole-genome Venn.  They ensure the factor pie totals and Venn
# circle totals describe exactly the same peak universe.
PEAK_LOCI: dict[str, Path] | None = None
BEDTOOLS = Path("bedtools")


def load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot import {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def merge_intervals(rows: pd.DataFrame) -> pd.DataFrame:
    """Match the r23 IgG display input: the union of its two strict calls."""
    merged: list[tuple[str, int, int]] = []
    for chrom, group in rows.sort_values(["chrom", "start", "end"]).groupby("chrom", sort=False):
        start: int | None = None
        end: int | None = None
        for row in group.itertuples(index=False):
            if start is None or int(row.start) > end:
                if start is not None:
                    merged.append((str(chrom), start, end))
                start, end = int(row.start), int(row.end)
            else:
                end = max(end, int(row.end))
        if start is not None:
            merged.append((str(chrom), start, end))
    return pd.DataFrame(merged, columns=["chrom", "start", "end"])


def materialize_igg_union() -> Path:
    if IGG_BED is not None:
        if not IGG_BED.is_file() or IGG_BED.stat().st_size == 0:
            raise FileNotFoundError(IGG_BED)
        table = pd.read_csv(IGG_BED, sep="\t", header=None, usecols=[0, 1, 2], names=["chrom", "start", "end"])
        table["start"] = pd.to_numeric(table["start"], errors="raise").astype(int)
        table["end"] = pd.to_numeric(table["end"], errors="raise").astype(int)
        table = table.loc[table["end"] > table["start"]].drop_duplicates()
        union = merge_intervals(table)
        output = STAGE / "IgG_2020_union_strict_peaks.bed"
        output.parent.mkdir(parents=True, exist_ok=True)
        union.to_csv(output, sep="\t", header=False, index=False)
        return output
    inputs = []
    for path in IGG_CALLS:
        if not path.is_file() or path.stat().st_size == 0:
            raise FileNotFoundError(path)
        inputs.append(pd.read_csv(path, sep="\t", header=None, usecols=[0, 1, 2], names=["chrom", "start", "end"]))
    table = pd.concat(inputs, ignore_index=True)
    table["start"] = pd.to_numeric(table["start"], errors="raise").astype(int)
    table["end"] = pd.to_numeric(table["end"], errors="raise").astype(int)
    table = table.loc[table["end"] > table["start"]].drop_duplicates()
    union = merge_intervals(table)
    output = STAGE / "IgG_2020_union_strict_peaks.bed"
    output.parent.mkdir(parents=True, exist_ok=True)
    union.to_csv(output, sep="\t", header=False, index=False)
    return output


def materialize_headerless_peak_bed(source: Path, label: str) -> Path:
    """Remove the tabular header from r23 interval exports before annotation."""
    table = pd.read_csv(source, sep="\t", header=None, names=["chrom", "start", "end"], dtype=str)
    table["start"] = pd.to_numeric(table["start"], errors="coerce")
    table["end"] = pd.to_numeric(table["end"], errors="coerce")
    table = table.dropna(subset=["start", "end"]).copy()
    table["start"] = table["start"].astype(int)
    table["end"] = table["end"].astype(int)
    table = table.loc[table["end"] > table["start"], ["chrom", "start", "end"]]
    output = STAGE / f"{label}.headerless.bed"
    table.to_csv(output, sep="\t", header=False, index=False)
    return output


def configure_renderer(peak_beds: dict[str, Path], unit_label: str):
    renderer = load_module(PROJECT / "peak_annotation.py", "peak_annotation")
    renderer.WORK_ROOT = CUTRUN
    renderer.OUTDIR = STAGE
    renderer.TMPDIR = STAGE / "intermediate"
    renderer.PAPER_TMPDIR = STAGE / "paperfig_intermediate"
    renderer.GTF_CANDIDATES = [GTF]
    renderer.ASSAYS = list(peak_beds)
    renderer.PEAK_UNIT_LABEL = unit_label
    renderer.USE_NT_KD_UNION = False
    renderer.BEDTOOLS = BEDTOOLS
    renderer.find_peaks_bed = lambda assay: peak_beds[assay]
    renderer.ensure_dirs()
    return renderer


def category_beds(renderer, first_peak: Path) -> dict[str, Path]:
    gtf = renderer.autodetect_gtf()
    annotation = renderer.harmonize_chrom_style_from_bed(first_peak, renderer.gtf_parse(gtf))
    transcripts = renderer.build_transcript_tss_tes(annotation, renderer.pick_canonical_transcripts(annotation))
    exons = renderer.build_transcript_exons(annotation, renderer.pick_canonical_transcripts(annotation))
    utr5, utr3 = renderer.build_transcript_utr(annotation, renderer.pick_canonical_transcripts(annotation), transcripts)
    if transcripts.empty:
        raise RuntimeError("No canonical transcripts were recovered from GENCODE M25")

    categories: dict[str, Path] = {}
    beds = {
        "Promoter (±1kb)": renderer.build_promoter_bed(transcripts),
        "5' UTR": utr5,
        "3' UTR": utr3,
        "Downstream (<=300bp)": renderer.build_downstream_bed(transcripts),
    }
    first_exon, other_exon = renderer.build_first_other_exons(exons)
    first_intron, other_intron = renderer.build_introns_from_exons(exons)
    beds.update({
        "1st Exon": first_exon,
        "Other Exon": other_exon,
        "1st Intron": first_intron,
        "Other Intron": other_intron,
    })
    for name, table in beds.items():
        path = renderer.TMPDIR / f"category_{name.replace(' ', '_').replace('/', '_')}.bed"
        renderer.write_bed3(table, path)
        categories[name] = path
    distal = renderer.TMPDIR / "category_Distal_Intergenic_placeholder.bed"
    distal.write_text("")
    categories["Distal Intergenic"] = distal
    return categories


def annotate_and_plot(peak_beds: dict[str, Path], unit_label: str, combined: bool) -> pd.DataFrame:
    renderer = configure_renderer(peak_beds, unit_label)
    categories = category_beds(renderer, next(iter(peak_beds.values())))
    results: dict[str, dict[str, int]] = {}
    totals: dict[str, int] = {}
    rows: list[dict[str, object]] = []
    for factor, peak_bed in peak_beds.items():
        midpoint_bed = renderer.TMPDIR / f"{factor}.peaks_mid.bed"
        total = renderer.bed_to_peak_midpoints(peak_bed, midpoint_bed)
        counts = renderer.assign_categories(midpoint_bed, categories) if total else {name: 0 for name in renderer.CATEGORY_ORDER}
        renderer.save_single(STAGE / f"FigX_peak_distribution_{factor}.png", factor, counts, total)
        results[factor] = counts
        totals[factor] = total
        rows.extend({"factor": factor, "category": category, "n_peaks": int(value), "total_peaks": total} for category, value in counts.items())
    if combined:
        renderer.save_three_panel(STAGE / "FigX_peak_distribution_pies.png", results, totals)
    return pd.DataFrame(rows)


def move_visuals(names: tuple[str, ...]) -> None:
    VISUALS.mkdir(parents=True, exist_ok=True)
    for name in names:
        source = STAGE / name
        if not source.is_file() or source.stat().st_size == 0:
            raise FileNotFoundError(source)
        source.replace(VISUALS / name)


def main() -> None:
    STAGE.mkdir(parents=True, exist_ok=True)
    factor_sources = PEAK_LOCI or {
        factor: PEAKS / f"{factor}_pooled_six_library_support4of6.bed"
        for factor in FACTORS
    }
    igg_inputs = (IGG_BED,) if IGG_BED is not None else IGG_CALLS
    missing = [str(path) for path in (*factor_sources.values(), *igg_inputs, GTF) if path is None or not path.is_file() or path.stat().st_size == 0]
    if missing:
        raise FileNotFoundError("Missing peak-distribution input(s):\n" + "\n".join(missing))
    factor_beds = {
        factor: materialize_headerless_peak_bed(source, factor)
        for factor, source in factor_sources.items()
    }
    igg_bed = materialize_igg_union()

    factor_counts = annotate_and_plot(factor_beds, "peaks", combined=True)
    igg_counts = annotate_and_plot({"IgG": igg_bed}, "peaks", combined=False)
    move_visuals((
        "FigX_peak_distribution_pies.png",
        "FigX_peak_distribution_MCM3.png",
        "FigX_peak_distribution_NONO.png",
        "FigX_peak_distribution_PSPC1.png",
        "FigX_peak_distribution_IgG.png",
    ))
    summary = pd.concat((factor_counts, igg_counts), ignore_index=True)
    summary.to_csv(DATA / "peak_distribution" / "Pie_PeakDistribution_counts.tsv", sep="\t", index=False)
    print(f"[DONE] Peak-distribution pies: {VISUALS}")


if __name__ == "__main__":
    main()
