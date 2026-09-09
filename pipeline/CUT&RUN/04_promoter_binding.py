#!/usr/bin/env python3
"""Call strict peaks and derive final 4-of-6 promoter-bound gene sets."""

from __future__ import annotations

import argparse
import math
import os
import shutil
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

import pandas as pd

PIPELINE_ROOT = Path(__file__).resolve().parents[1]
if str(PIPELINE_ROOT) not in sys.path:
    sys.path.insert(0, str(PIPELINE_ROOT))

from config import (
    CUTRUN_ROOT, FACTORS, PEAK_FE_MIN, PEAK_P_MAX, PEAK_Q_2019_MAX,
    PEAK_Q_2020_MAX, PROMOTER_BP, PROMOTER_OVERLAP_BP, RATIO_MIN,
    REFERENCE_DIR, TOTAL_FACTOR_SUPPORT,
)


def executable(name: str, candidates: tuple[str, ...]) -> str:
    for candidate in candidates:
        hit = shutil.which(candidate)
        if hit:
            return hit
        if Path(candidate).is_file():
            return candidate
    raise FileNotFoundError(f"Required executable not found: {name}")


def call_macs(macs3: str, treatment: list[Path], output_dir: Path, name: str, control: list[Path] | None, dry_run: bool) -> Path:
    output_dir.mkdir(parents=True, exist_ok=True)
    peak = output_dir / f"{name}_peaks.narrowPeak"
    if peak.is_file():
        return peak
    command = [macs3, "callpeak", "-t", *(str(path) for path in treatment)]
    if control:
        command += ["-c", *(str(path) for path in control)]
    command += ["-f", "BAMPE", "-g", "mm", "--keep-dup", "all", "--call-summits", "-p", f"{PEAK_P_MAX:g}", "-n", name, "--outdir", str(output_dir)]
    print("[RUN]", " ".join(command), flush=True)
    if not dry_run:
        subprocess.run(command, check=True)
    return peak


def read_filtered_peaks(path: Path, q_max: float) -> list[tuple[str, int, int]]:
    if not path.is_file():
        raise FileNotFoundError(path)
    p_min = -math.log10(PEAK_P_MAX)
    q_min = -math.log10(q_max)
    result: list[tuple[str, int, int]] = []
    with path.open() as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            chrom, start, end = fields[0], int(fields[1]), int(fields[2])
            fe, minus_log10_p, minus_log10_q = map(float, fields[6:9])
            if end > start and fe >= PEAK_FE_MIN and minus_log10_p >= p_min and minus_log10_q >= q_min:
                result.append((chrom, start, end))
    return result


def merge_intervals(rows: list[tuple[str, int, int]]) -> list[tuple[str, int, int]]:
    merged: list[tuple[str, int, int]] = []
    for chrom, start, end in sorted(rows):
        if not merged or chrom != merged[-1][0] or start > merged[-1][2]:
            merged.append((chrom, start, end))
        else:
            old_chrom, old_start, old_end = merged[-1]
            merged[-1] = old_chrom, old_start, max(old_end, end)
    return merged


def supported_segments(calls: list[list[tuple[str, int, int]]], minimum_support: int) -> list[tuple[str, int, int]]:
    events: dict[str, list[tuple[int, int, int]]] = defaultdict(list)
    for index, call in enumerate(calls):
        bit = 1 << index
        for chrom, start, end in call:
            events[chrom].append((start, 1, bit)); events[chrom].append((end, -1, bit))
    result: list[tuple[str, int, int]] = []
    for chrom, chrom_events in events.items():
        active = 0; previous: int | None = None
        for position, kind, bit in sorted(chrom_events, key=lambda row: (row[0], row[1])):
            if previous is not None and position > previous and active.bit_count() >= minimum_support:
                result.append((chrom, previous, position))
            active = (active | bit) if kind == 1 else (active & ~bit)
            previous = position
    return merge_intervals(result)


def gtf_promoters(path: Path) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    with path.open() as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            field = line.rstrip("\n").split("\t")
            if len(field) != 9 or field[2] != "gene":
                continue
            attrs = field[8]
            gene = next((item.split('"')[1] for item in attrs.split("; ") if item.startswith("gene_name ")), None)
            if not gene:
                continue
            start, end, strand = int(field[3]) - 1, int(field[4]), field[6]
            tss = start if strand == "+" else end - 1
            rows.append({"chrom": field[0], "start": max(0, tss - PROMOTER_BP), "end": tss + PROMOTER_BP + 1, "gene": gene, "strand": strand})
    table = pd.DataFrame(rows).drop_duplicates(["chrom", "start", "end", "gene"])
    if table.empty:
        raise RuntimeError("No GENCODE M25 promoters parsed")
    return table


def overlap_genes(peaks: list[tuple[str, int, int]], promoters: pd.DataFrame) -> set[str]:
    by_chrom: dict[str, list[tuple[int, int, str]]] = defaultdict(list)
    for row in promoters.itertuples(index=False):
        by_chrom[str(row.chrom)].append((int(row.start), int(row.end), str(row.gene)))
    hits: set[str] = set()
    for chrom, start, end in peaks:
        for p_start, p_end, gene in by_chrom.get(chrom, []):
            if min(end, p_end) - max(start, p_start) >= PROMOTER_OVERLAP_BP:
                hits.add(gene)
    return hits


def represented_mean(bigwig: Path, chrom: str, start: int, end: int) -> float:
    import pyBigWig
    with pyBigWig.open(str(bigwig)) as reader:
        intervals = reader.intervals(chrom, start, end) or []
    numerator = 0.0; denominator = 0
    for left, right, value in intervals:
        overlap = max(0, min(end, int(right)) - max(start, int(left)))
        if overlap:
            numerator += overlap * float(value); denominator += overlap
    return numerator / denominator if denominator else 0.0


def promoter_ratio_support(promoters: pd.DataFrame, factor_tracks: list[Path], igg_tracks: list[Path]) -> pd.DataFrame:
    rows = []
    for row in promoters.itertuples(index=False):
        ratios = []
        for factor_bw, igg_bw in zip(factor_tracks, igg_tracks):
            factor = represented_mean(factor_bw, row.chrom, int(row.start), int(row.end))
            igg = represented_mean(igg_bw, row.chrom, int(row.start), int(row.end))
            ratio = float("inf") if igg == 0 and factor > 0 else (factor / igg if igg > 0 else 0.0)
            ratios.append(ratio)
        rows.append({"gene": row.gene, "ratio_rep1": ratios[0], "ratio_rep2": ratios[1], "ratio_support": ratios[0] >= RATIO_MIN and ratios[1] >= RATIO_MIN})
    return pd.DataFrame(rows)


def write_bed(rows: list[tuple[str, int, int]], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(path, sep="\t", header=False, index=False)


def copy_file(source: Path, destination: Path) -> None:
    if not source.is_file():
        raise FileNotFoundError(f"Missing accepted membership input: {source}")
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, destination)


def write_presentation_gene_beds() -> None:
    """Create GENCODE M25 gene-body and promoter coordinate files."""
    gtf = REFERENCE_DIR / "gencode.vM25.annotation.gtf"
    if not gtf.is_file():
        raise FileNotFoundError(f"Missing GENCODE M25 annotation: {gtf}")
    bodies: list[tuple[str, int, int, str, int, str]] = []
    promoters: list[tuple[str, int, int, str]] = []
    with gtf.open() as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] != "gene":
                continue
            name = next((item.split('"')[1] for item in fields[8].split("; ") if item.startswith("gene_name ")), None)
            if not name:
                continue
            start, end, strand = int(fields[3]) - 1, int(fields[4]), fields[6]
            tss = start if strand == "+" else end - 1
            bodies.append((fields[0], start, end, name, 0, strand))
            promoters.append((fields[0], max(0, tss - PROMOTER_BP), tss + PROMOTER_BP + 1, name))
    data = CUTRUN_ROOT / "data"
    data.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(bodies).drop_duplicates(3).to_csv(data / "GeneBodies_M25.bed6", sep="\t", header=False, index=False)
    pd.DataFrame(promoters).drop_duplicates(3).to_csv(data / "Promoters_M25_TSSplusminus1kb.bed", sep="\t", header=False, index=False)


REGION_NAMES = {
    1: "100_only_MCM3", 2: "010_only_NONO", 4: "001_only_PSPC1",
    3: "110_MCM3_NONO", 5: "101_MCM3_PSPC1", 6: "011_NONO_PSPC1",
    7: "111_MCM3_NONO_PSPC1",
}


def peak_locus_membership(sets: dict[str, list[tuple[str, int, int]]]) -> pd.DataFrame:
    """Return the canonical nonredundant locus universe for all peak panels.

    Peak calls from the three factors are merged into connected genomic loci.
    A locus is present in a factor when one or more retained consensus peaks
    from that factor overlap the locus.  The returned table is the one source
    of membership truth for the whole-genome Venn and all three factor peak
    distribution pies.  Promoter membership remains independently defined by
    the stricter >=250-bp peak-promoter rule.
    """
    all_intervals = sorted(
        interval
        for factor in FACTORS
        for interval in sets[factor]
    )
    loci: list[tuple[str, int, int]] = []
    for chrom, start, end in all_intervals:
        if loci and chrom == loci[-1][0] and start <= loci[-1][2]:
            loci[-1] = chrom, loci[-1][1], max(end, loci[-1][2])
        else:
            loci.append((chrom, start, end))

    by_chrom: dict[str, dict[str, list[tuple[int, int]]]] = {
        factor: defaultdict(list) for factor in FACTORS
    }
    for factor in FACTORS:
        for chrom, start, end in sets[factor]:
            by_chrom[factor][chrom].append((start, end))

    rows: list[dict[str, object]] = []
    for chrom, start, end in loci:
        mask = 0
        for index, factor in enumerate(FACTORS):
            if any(left < end and right > start for left, right in by_chrom[factor][chrom]):
                mask |= 1 << index
        if mask:
            rows.append({"region_mask": mask, "chrom": chrom, "start": start, "end": end})

    table = pd.DataFrame(rows, columns=["region_mask", "chrom", "start", "end"])
    if table.empty:
        return table.assign(locus_id=pd.Series(dtype=str), membership=pd.Series(dtype=str))
    table.insert(0, "locus_id", [f"PeakLocus_{index:06d}" for index in range(1, len(table) + 1)])
    table["membership"] = table["region_mask"].map(REGION_NAMES)
    return table


def peak_locus_region_counts(sets: dict[str, list[tuple[str, int, int]]]) -> pd.DataFrame:
    """Summarize the canonical peak-locus universe as Venn regions."""
    membership = peak_locus_membership(sets)
    rows = membership.loc[:, ["region_mask"]] if not membership.empty else pd.DataFrame(columns=["region_mask"])

    counts = pd.DataFrame(rows).groupby("region_mask", as_index=False).size().rename(columns={"size": "peak_loci"})
    counts = pd.DataFrame({"region_mask": range(1, 8)}).merge(counts, on="region_mask", how="left").fillna({"peak_loci": 0})
    counts["peak_loci"] = counts["peak_loci"].astype(int)
    counts["minimum_peak_overlap_bp"] = 1
    counts["unit"] = "merged whole-genome consensus peak locus"
    return counts


def write_peak_locus_inputs(
    sets: dict[str, list[tuple[str, int, int]]], output_dir: Path
) -> pd.DataFrame:
    """Materialize factor-associated loci used identically by Venn and pies."""
    output_dir.mkdir(parents=True, exist_ok=True)
    membership = peak_locus_membership(sets)
    membership.to_csv(output_dir / "Venn_Peaks_loci.tsv", sep="\t", index=False)
    for index, factor in enumerate(FACTORS):
        bit = 1 << index
        subset = membership.loc[
            membership["region_mask"].astype(int).map(lambda value: bool(value & bit)),
            ["chrom", "start", "end"],
        ]
        subset.to_csv(output_dir / f"Venn_Peaks_{factor}.bed", sep="\t", header=False, index=False)
    return membership


def write_peak_universe_summary(
    peak_sets: dict[str, list[tuple[str, int, int]]],
    promoter_sets: dict[str, set[str]],
    output: Path,
    igg_peaks: list[tuple[str, int, int]] | None = None,
) -> None:
    """Record the distinct units behind peak, promoter, and profile panels."""
    membership = peak_locus_membership(peak_sets)
    rows: list[dict[str, object]] = []
    for index, factor in enumerate(FACTORS):
        bit = 1 << index
        n_loci = int(membership["region_mask"].astype(int).map(lambda value: bool(value & bit)).sum())
        rows.extend((
            {"panel": "Venn_Peaks", "set": factor, "unit": "factor-associated merged consensus peak loci", "n": n_loci},
            {"panel": "Pie_PeakDistribution", "set": factor, "unit": "factor-associated merged consensus peak loci", "n": n_loci},
            {"panel": "Venn_PromoterGenes", "set": factor, "unit": "promoter-bound genes", "n": len(promoter_sets[factor])},
        ))
    rows.extend((
        {"panel": "Profile_MCM3_NONO_PSPC1", "set": "MCM3+NONO+PSPC1", "unit": "co-bound promoter genes", "n": len(set.intersection(*promoter_sets.values()))},
        {"panel": "Profile_MCM3_NONO", "set": "MCM3+NONO only", "unit": "co-bound promoter genes", "n": len((promoter_sets["MCM3"] & promoter_sets["NONO"]) - promoter_sets["PSPC1"])},
        {"panel": "Profile_MCM3_PSPC1", "set": "MCM3+PSPC1 only", "unit": "co-bound promoter genes", "n": len((promoter_sets["MCM3"] & promoter_sets["PSPC1"]) - promoter_sets["NONO"])},
    ))
    if igg_peaks is not None:
        rows.append({
            "panel": "Pie_IgG",
            "set": "IgG",
            "unit": "merged IgG control peak loci",
            "n": len(merge_intervals(igg_peaks)),
        })
    pd.DataFrame(rows).to_csv(output, sep="\t", index=False)


def read_bed_intervals(path: Path) -> list[tuple[str, int, int]]:
    """Read BED-like intervals while tolerating a single header row."""
    intervals: list[tuple[str, int, int]] = []
    with path.open() as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                continue
            try:
                start, end = int(fields[1]), int(fields[2])
            except ValueError:
                continue
            if end > start:
                intervals.append((fields[0], start, end))
    return intervals


def write_raw_renderer_inputs(
    final_sets: dict[str, set[str]],
    final_peaks: dict[str, list[tuple[str, int, int]]],
    peak_summary: list[dict[str, object]],
) -> None:
    """Build the source-layout tables consumed by the locked figure renderers.

    This compatibility layout is generated from the fresh raw peak and promoter
    products.  It is deliberately equivalent to the locally materialized
    accepted-input layout used by the fast path, not an import from another
    reanalysis directory.
    """
    source = CUTRUN_ROOT / "data" / "figure_inputs"
    promoter = source / "promoter_gene_venn"
    peaks = source / "batch_stratified_peaks"
    promoter.mkdir(parents=True, exist_ok=True)
    peaks.mkdir(parents=True, exist_ok=True)
    (source / "peak_distribution").mkdir(parents=True, exist_ok=True)

    for factor in FACTORS:
        genes = pd.DataFrame({"gene": sorted(final_sets[factor])})
        genes.to_csv(promoter / f"Venn_PromoterGenes_{factor}.tsv", sep="\t", index=False)
        support = CUTRUN_ROOT / "05_promoters" / f"Venn_PromoterGenes_{factor}_support.tsv"
        if not support.is_file():
            raise FileNotFoundError(support)
        shutil.copy2(support, promoter / f"Venn_PromoterGenes_{factor}_support.tsv")
        write_bed(final_peaks[factor], peaks / f"Venn_Peaks_{factor}_support4of6.bed")

    regions: list[dict[str, str]] = []
    for gene in sorted(set().union(*final_sets.values())):
        mask = sum((1 << index) for index, factor in enumerate(FACTORS) if gene in final_sets[factor])
        regions.append({"region": REGION_NAMES[mask], "gene": gene})
    pd.DataFrame(regions).to_csv(
        promoter / "Venn_PromoterGenes_regions.tsv", sep="\t", index=False
    )
    parameters = CUTRUN_ROOT / "data" / "AnalysisParameters.tsv"
    if parameters.is_file():
        shutil.copy2(parameters, promoter / "AnalysisParameters.tsv")

    peak_locus_region_counts(final_peaks).to_csv(
        source / "Venn_Peaks_counts.tsv", sep="\t", index=False
    )
    write_peak_locus_inputs(final_peaks, source / "peak_loci")
    write_peak_universe_summary(final_peaks, final_sets, source / "Venn_Peaks_universe.tsv")
    summary = pd.DataFrame(peak_summary).rename(columns={"peak_segments": "n_pooled_six_library_support_segments"})
    summary.to_csv(source / "Venn_Peaks_support.tsv", sep="\t", index=False)

    igg_source = CUTRUN_ROOT / "data" / "PeakSets" / "Peaks_IgG.bed"
    if not igg_source.is_file():
        raise FileNotFoundError(igg_source)
    shutil.copy2(igg_source, source / "peak_distribution" / "Peaks_IgG.bed")


def materialize_accepted_membership() -> None:
    """Stage locally copied 4-of-6 peak and promoter evidence for rendering.

    The accepted membership inputs were copied by
    ``01_prepare.py``. This stage converts their long,
    analysis-specific filenames into concise final data names.  It never reads
    a previous reanalysis directory.
    """
    source = CUTRUN_ROOT / "data" / "figure_inputs"
    promoter_source = source / "promoter_gene_venn"
    peak_source = source / "batch_stratified_peaks"
    data = CUTRUN_ROOT / "data"
    promoter_dir = CUTRUN_ROOT / "05_promoters"
    peak_dir = data / "PeakSets"
    data.mkdir(parents=True, exist_ok=True)
    promoter_dir.mkdir(parents=True, exist_ok=True)
    peak_dir.mkdir(parents=True, exist_ok=True)
    summary_rows: list[dict[str, object]] = []
    for factor in FACTORS:
        genes = promoter_source / f"Venn_PromoterGenes_{factor}.tsv"
        support = promoter_source / f"Venn_PromoterGenes_{factor}_support.tsv"
        bed = peak_source / f"Venn_Peaks_{factor}_support4of6.bed"
        copy_file(genes, promoter_dir / f"Venn_PromoterGenes_{factor}.tsv")
        copy_file(support, promoter_dir / f"Venn_PromoterGenes_{factor}_support.tsv")
        copy_file(bed, peak_dir / f"Peaks_{factor}.bed")
        summary_rows.append({
            "factor": factor,
            "promoter_bound_genes": len(pd.read_csv(genes, sep="\t")),
            "peak_segments": sum(1 for _ in bed.open()),
            "rule": "total support >=4 of 6 independent factor libraries",
        })
    # This strict IgG union is used only for the IgG peak-annotation pie.
    copy_file(
        source / "peak_distribution" / "Peaks_IgG.bed",
        peak_dir / "Peaks_IgG.bed",
    )
    peak_sets = {
        factor: read_bed_intervals(peak_dir / f"Peaks_{factor}.bed")
        for factor in FACTORS
    }
    peak_locus_region_counts(peak_sets).to_csv(
        data / "Venn_Peaks_counts.tsv", sep="\t", index=False
    )
    peak_locus_region_counts(peak_sets).to_csv(
        source / "Venn_Peaks_counts.tsv", sep="\t", index=False
    )
    promoter_sets = {
        factor: set(pd.read_csv(promoter_dir / f"Venn_PromoterGenes_{factor}.tsv", sep="\t", dtype=str)["gene"].dropna())
        for factor in FACTORS
    }
    write_peak_locus_inputs(peak_sets, data / "PeakLoci")
    write_peak_locus_inputs(peak_sets, source / "peak_loci")
    igg_peaks = read_bed_intervals(peak_dir / "Peaks_IgG.bed")
    write_peak_universe_summary(peak_sets, promoter_sets, data / "Venn_Peaks_universe.tsv", igg_peaks)
    copy_file(source / "Venn_Peaks_support.tsv", data / "Venn_Peaks_support.tsv")
    copy_file(promoter_source / "AnalysisParameters.tsv", data / "AnalysisParameters.tsv")
    pd.DataFrame(summary_rows).to_csv(data / "Venn_PromoterGenes_summary.tsv", sep="\t", index=False)
    track_manifest = CUTRUN_ROOT / "03_bigwig" / "IGV_representation_libraries.tsv"
    if track_manifest.is_file():
        table = pd.read_csv(track_manifest, sep="\t")
        table.to_csv(CUTRUN_ROOT / "03_bigwig" / "YeastNormalization.tsv", sep="\t", index=False)
    write_presentation_gene_beds()
    print("[DONE] Staged accepted six-library peak and promoter membership inputs locally.")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--from-raw", action="store_true", help="Rebuild peak and promoter membership from local raw-derived BAMs")
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()
    manifest_path = CUTRUN_ROOT / "metadata" / "Samples.tsv"
    if not manifest_path.is_file():
        raise FileNotFoundError("Run 01_prepare.py first")
    manifest = pd.read_csv(manifest_path, sep="\t")
    if not args.from_raw:
        if args.dry_run:
            print("[DRY-RUN] Would stage locally materialized accepted peak and promoter membership inputs.")
            return
        materialize_accepted_membership()
        return
    macs3 = executable("macs3", ("macs3", "/opt/anaconda3/envs/cutrun_env/bin/macs3"))
    if args.dry_run:
        print(
            "[DRY-RUN] Would call 2019 individual-factor MACS3 peaks and 2020 "
            "matched-IgG pooled/replicate peaks, then apply the configured promoter rule."
        )
        return
    promoters = gtf_promoters(REFERENCE_DIR / "gencode.vM25.annotation.gtf")
    data = CUTRUN_ROOT / "data"; peak_root = CUTRUN_ROOT / "04_peaks"; promoter_dir = CUTRUN_ROOT / "05_promoters"
    data.mkdir(parents=True, exist_ok=True); promoter_dir.mkdir(parents=True, exist_ok=True)
    write_presentation_gene_beds()

    all_support_rows: list[dict[str, object]] = []
    final_sets: dict[str, set[str]] = {}
    final_peaks: dict[str, list[tuple[str, int, int]]] = {}
    for factor in FACTORS:
        rows19 = manifest.loc[(manifest.batch == "2019") & (manifest.factor == factor) & manifest.selected_factor_library]
        rows20 = manifest.loc[(manifest.batch == "2020") & (manifest.factor == factor)]
        igg20 = manifest.loc[(manifest.batch == "2020") & (manifest.factor == "IgG")].sort_values("replicate")
        if len(rows19) != 4 or len(rows20) != 2 or len(igg20) != 2:
            raise RuntimeError(f"Unexpected library count for {factor}")
        calls19 = []
        support19: dict[str, int] = defaultdict(int)
        for row in rows19.itertuples(index=False):
            bam = CUTRUN_ROOT / "02_bam" / f"{row.library_id}.bam"
            raw = call_macs(macs3, [bam], peak_root / "2019" / factor / str(row.library_id), str(row.library_id), None, args.dry_run)
            if args.dry_run:
                continue
            passed = read_filtered_peaks(raw, PEAK_Q_2019_MAX); calls19.append(passed)
            for gene in overlap_genes(passed, promoters):
                support19[gene] += 1
        if args.dry_run:
            continue
        factor_bams20 = [CUTRUN_ROOT / "02_bam" / f"{row.library_id}.bam" for row in rows20.sort_values("replicate").itertuples(index=False)]
        igg_bams20 = [CUTRUN_ROOT / "02_bam" / f"{row.library_id}.bam" for row in igg20.itertuples(index=False)]
        raw_pooled = call_macs(macs3, factor_bams20, peak_root / "2020" / factor / "pooled", f"{factor}_pooled", igg_bams20, args.dry_run)
        pooled = read_filtered_peaks(raw_pooled, PEAK_Q_2020_MAX)
        candidate20 = overlap_genes(pooled, promoters)
        factor_tracks20 = [CUTRUN_ROOT / "03_bigwig" / str(row.library_id) / f"{row.library_id}.bw" for row in rows20.sort_values("replicate").itertuples(index=False)]
        igg_tracks20 = [CUTRUN_ROOT / "03_bigwig" / str(row.library_id) / f"{row.library_id}.bw" for row in igg20.itertuples(index=False)]
        ratios = promoter_ratio_support(promoters.loc[promoters.gene.isin(candidate20)], factor_tracks20, igg_tracks20)
        support20 = set(ratios.loc[ratios.ratio_support, "gene"])
        final = {gene for gene in set(support19) | candidate20 if support19.get(gene, 0) + (2 if gene in support20 else 0) >= TOTAL_FACTOR_SUPPORT}
        final_sets[factor] = final
        support_table = pd.DataFrame({"gene": sorted(set(support19) | candidate20)})
        support_table["n_2019_libraries"] = support_table.gene.map(support19).fillna(0).astype(int)
        support_table = support_table.merge(ratios, on="gene", how="left").fillna({"ratio_support": False})
        support_table["n_2020_libraries"] = support_table.ratio_support.astype(int) * 2
        support_table["n_total_libraries"] = support_table.n_2019_libraries + support_table.n_2020_libraries
        support_table["promoter_bound"] = support_table.n_total_libraries >= TOTAL_FACTOR_SUPPORT
        support_table.to_csv(promoter_dir / f"Venn_PromoterGenes_{factor}_support.tsv", sep="\t", index=False)
        pd.DataFrame({"gene": sorted(final)}).to_csv(promoter_dir / f"Venn_PromoterGenes_{factor}.tsv", sep="\t", index=False)
        # Whole-genome figure segments are the same 4-of-6 atomic support rule.
        calls20 = []
        for row, control in zip(rows20.sort_values("replicate").itertuples(index=False), igg20.itertuples(index=False)):
            raw = call_macs(macs3, [CUTRUN_ROOT / "02_bam" / f"{row.library_id}.bam"], peak_root / "2020" / factor / str(row.library_id), str(row.library_id), [CUTRUN_ROOT / "02_bam" / f"{control.library_id}.bam"], args.dry_run)
            calls20.append(read_filtered_peaks(raw, PEAK_Q_2020_MAX))
        final_peaks[factor] = supported_segments(calls19 + calls20, TOTAL_FACTOR_SUPPORT)
        write_bed(final_peaks[factor], data / "PeakSets" / f"Peaks_{factor}.bed")
        all_support_rows.append({
            "factor": factor,
            "n_2019_calls": sum(len(call) for call in calls19),
            "n_2020_calls": len(pooled),
            "promoter_bound_genes": len(final),
            "peak_segments": len(final_peaks[factor]),
            "minimum_total_support": TOTAL_FACTOR_SUPPORT,
            "rule": "4 of 6 independent factor libraries",
        })
    if args.dry_run:
        return
    # Strict IgG calls are retained for the IgG distribution pie only.
    igg_rows = manifest.loc[(manifest.batch == "2020") & (manifest.factor == "IgG")]
    igg_calls = []
    for row in igg_rows.itertuples(index=False):
        raw = call_macs(macs3, [CUTRUN_ROOT / "02_bam" / f"{row.library_id}.bam"], peak_root / "2020" / "IgG" / str(row.library_id), str(row.library_id), None, False)
        igg_calls.append(read_filtered_peaks(raw, PEAK_Q_2020_MAX))
    igg_peaks = merge_intervals([peak for call in igg_calls for peak in call])
    write_bed(igg_peaks, data / "PeakSets" / "Peaks_IgG.bed")
    pd.DataFrame(all_support_rows).to_csv(data / "Venn_PromoterGenes_summary.tsv", sep="\t", index=False)
    pd.DataFrame([{"parameter": "promoter", "value": "GENCODE M25 TSS +/-1,000 bp"}, {"parameter": "overlap", "value": ">=250 bp"}, {"parameter": "final_support", "value": ">=4/6 independent factor libraries"}, {"parameter": "2019_peak", "value": "p<=1e-4; q<=0.05; FE>=3"}, {"parameter": "2020_peak", "value": "pooled matched IgG p<=1e-4; q<=0.01; FE>=3; ratio>=2 in both replicates"}]).to_csv(data / "AnalysisParameters.tsv", sep="\t", index=False)
    peak_locus_region_counts(final_peaks).to_csv(data / "Venn_Peaks_counts.tsv", sep="\t", index=False)
    write_peak_locus_inputs(final_peaks, data / "PeakLoci")
    write_peak_universe_summary(final_peaks, final_sets, data / "Venn_Peaks_universe.tsv", igg_peaks)
    write_raw_renderer_inputs(final_sets, final_peaks, all_support_rows)


if __name__ == "__main__":
    main()
