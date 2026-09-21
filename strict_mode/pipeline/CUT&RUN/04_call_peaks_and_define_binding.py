#!/usr/bin/env python3
"""Call candidate peaks, apply final filters, and define promoter binding."""

from __future__ import annotations

import argparse
import json
import math
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
    CUTRUN_ROOT,
    FACTORS,
    PEAK_CANDIDATE_P_MAX,
    PEAK_FE_MIN,
    PEAK_P_MAX,
    PEAK_Q_2020_MAX,
    PROMOTER_BP,
    PROMOTER_OVERLAP_BP,
    RATIO_MIN,
    REFERENCE_DIR,
)

REGION_NAMES = {
    1: "100_only_MCM3",
    2: "010_only_NONO",
    4: "001_only_PSPC1",
    3: "110_MCM3_NONO",
    5: "101_MCM3_PSPC1",
    6: "011_NONO_PSPC1",
    7: "111_MCM3_NONO_PSPC1",
}


def executable(name: str, candidates: tuple[str, ...]) -> str:
    for candidate in candidates:
        hit = shutil.which(candidate)
        if hit:
            return hit
        if Path(candidate).is_file():
            return candidate
    raise FileNotFoundError(f"Required executable not found: {name}")


def call_macs(
    macs3: str,
    treatment: list[Path],
    output: Path,
    name: str,
    controls: list[Path] | None,
) -> Path:
    output.mkdir(parents=True, exist_ok=True)
    peak = output / f"{name}_peaks.narrowPeak"
    command = [macs3, "callpeak", "-t", *(str(path) for path in treatment)]
    if controls:
        command.extend(("-c", *(str(path) for path in controls)))
    command.extend(
        (
            "-f",
            "BAMPE",
            "-g",
            "mm",
            "--keep-dup",
            "all",
            "--call-summits",
            "-p",
            f"{PEAK_CANDIDATE_P_MAX:g}",
            "-n",
            name,
            "--outdir",
            str(output),
        )
    )
    inputs = [*treatment, *(controls or [])]
    provenance = {
        "command": command,
        "inputs": [
            {
                "path": str(path.resolve()),
                "size": path.stat().st_size,
                "mtime_ns": path.stat().st_mtime_ns,
            }
            for path in inputs
        ],
    }
    provenance_path = output / f"{name}.provenance.json"
    if peak.is_file() and provenance_path.is_file():
        try:
            if json.loads(provenance_path.read_text()) == provenance:
                return peak
        except (json.JSONDecodeError, OSError):
            pass
    for suffix in ("_peaks.narrowPeak", "_peaks.xls", "_summits.bed", "_model.r"):
        (output / f"{name}{suffix}").unlink(missing_ok=True)
    print("[RUN]", " ".join(command), flush=True)
    subprocess.run(command, check=True)
    if not peak.is_file() or peak.stat().st_size == 0:
        raise RuntimeError(f"MACS3 did not produce a non-empty peak file: {peak}")
    provenance_path.write_text(json.dumps(provenance, indent=2) + "\n")
    return peak


def peak_passes(fields: list[str]) -> bool:
    p_min = -math.log10(PEAK_P_MAX)
    q_min = -math.log10(PEAK_Q_2020_MAX)
    if len(fields) < 9:
        return False
    start, end = int(fields[1]), int(fields[2])
    fold_enrichment, p_score, q_score = map(float, fields[6:9])
    return (
        end > start
        and p_score >= p_min
        and q_score >= q_min
        and fold_enrichment >= PEAK_FE_MIN
    )


def filter_peak_file(source: Path, destination: Path) -> Path:
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_name(destination.name + ".tmp")
    table = pd.read_csv(source, sep="\t", header=None)
    p_min = -math.log10(PEAK_P_MAX)
    q_min = -math.log10(PEAK_Q_2020_MAX)
    table = table[
        (pd.to_numeric(table.iloc[:, 7], errors="coerce") >= p_min)
        & (pd.to_numeric(table.iloc[:, 8], errors="coerce") >= q_min)
        & (pd.to_numeric(table.iloc[:, 6], errors="coerce") >= PEAK_FE_MIN)
    ]
    if table.empty:
        raise RuntimeError(f"Final peak filters removed every candidate peak: {source}")
    table.to_csv(temporary, sep="\t", header=False, index=False)
    temporary.replace(destination)
    return destination


def read_filtered_peaks(path: Path) -> list[tuple[str, int, int]]:
    peaks = []
    with path.open() as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if peak_passes(fields):
                peaks.append((fields[0], int(fields[1]), int(fields[2])))
    return peaks


def merge_intervals(rows: list[tuple[str, int, int]]) -> list[tuple[str, int, int]]:
    merged: list[tuple[str, int, int]] = []
    for chrom, start, end in sorted(rows):
        if not merged or chrom != merged[-1][0] or start > merged[-1][2]:
            merged.append((chrom, start, end))
        else:
            merged[-1] = chrom, merged[-1][1], max(merged[-1][2], end)
    return merged


def supported_segments(
    calls: list[list[tuple[str, int, int]]],
    minimum: int,
) -> list[tuple[str, int, int]]:
    events: dict[str, list[tuple[int, int, int]]] = defaultdict(list)
    for index, call in enumerate(calls):
        bit = 1 << index
        for chrom, start, end in call:
            events[chrom].extend(((start, 1, bit), (end, -1, bit)))
    segments = []
    for chrom, chrom_events in events.items():
        active = 0
        previous = None
        for position, kind, bit in sorted(
            chrom_events, key=lambda row: (row[0], row[1])
        ):
            if (
                previous is not None
                and position > previous
                and active.bit_count() >= minimum
            ):
                segments.append((chrom, previous, position))
            active = active | bit if kind == 1 else active & ~bit
            previous = position
    return merge_intervals(segments)


def atomic_peak_segments(
    peak_sets: dict[str, list[tuple[str, int, int]]]
) -> pd.DataFrame:
    events: dict[str, list[tuple[int, int, int]]] = defaultdict(list)
    for index, factor in enumerate(FACTORS):
        bit = 1 << index
        for chrom, start, end in peak_sets[factor]:
            events[chrom].extend(((start, 1, bit), (end, -1, bit)))
    rows = []
    for chrom, chrom_events in events.items():
        active = 0
        previous = None
        for position, kind, bit in sorted(
            chrom_events, key=lambda row: (row[0], row[1])
        ):
            if (
                previous is not None
                and position - previous >= PROMOTER_OVERLAP_BP
                and active
            ):
                rows.append((chrom, previous, position, active))
            active = active | bit if kind == 1 else active & ~bit
            previous = position
    table = pd.DataFrame(rows, columns=["chrom", "start", "end", "region_mask"])
    table = table.sort_values(["chrom", "start", "end"], kind="stable").reset_index(
        drop=True
    )
    table.insert(
        0,
        "locus_id",
        table.apply(
            lambda row: f"{row.chrom}:{int(row.start)}-{int(row.end)}", axis=1
        ),
    )
    table["membership"] = table.region_mask.map(REGION_NAMES)
    return table


def gtf_tables(path: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    promoters, bodies = [], []
    with path.open() as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] != "gene":
                continue
            gene = next(
                (
                    item.split('"')[1]
                    for item in fields[8].split("; ")
                    if item.startswith("gene_name ")
                ),
                None,
            )
            if not gene:
                continue
            start, end, strand = int(fields[3]) - 1, int(fields[4]), fields[6]
            tss = start if strand == "+" else end
            promoters.append(
                (fields[0], max(0, tss - PROMOTER_BP), tss + PROMOTER_BP, gene)
            )
            bodies.append((fields[0], start, end, gene, 0, strand))
    return (
        pd.DataFrame(
            promoters, columns=["chrom", "start", "end", "gene"]
        ).drop_duplicates(),
        pd.DataFrame(
            bodies, columns=["chrom", "start", "end", "gene", "score", "strand"]
        ).drop_duplicates("gene"),
    )


def overlap_genes(
    peaks: list[tuple[str, int, int]], promoters: pd.DataFrame
) -> set[str]:
    by_chrom = {
        chrom: list(
            group[["start", "end", "gene"]].itertuples(index=False, name=None)
        )
        for chrom, group in promoters.groupby("chrom", sort=False)
    }
    genes = set()
    for chrom, start, end in peaks:
        for p_start, p_end, gene in by_chrom.get(chrom, []):
            if min(end, p_end) - max(start, p_start) >= PROMOTER_OVERLAP_BP:
                genes.add(gene)
    return genes


def represented_mean(bigwig: Path, chrom: str, start: int, end: int) -> float:
    import pyBigWig

    with pyBigWig.open(str(bigwig)) as reader:
        if chrom not in reader.chroms() or start >= reader.chroms()[chrom]:
            return 0.0
        value = reader.stats(
            chrom,
            start,
            min(end, reader.chroms()[chrom]),
            type="mean",
            nBins=1,
            exact=True,
        )[0]
    return float(value or 0.0)


def promoter_support(
    promoters: pd.DataFrame,
    candidates: set[str],
    factor: str,
) -> pd.DataFrame:
    tracks = CUTRUN_ROOT / "03_bigwig"
    rows = []
    for gene, coordinates in promoters.loc[
        promoters.gene.isin(candidates)
    ].groupby("gene", sort=False):
        ratios = []
        for replicate in (1, 2):
            factor_bw = tracks / f"{factor}-{replicate}" / f"{factor}-{replicate}.bw"
            igg_bw = tracks / f"IgG-{replicate}" / f"IgG-{replicate}.bw"
            best = 0.0
            for row in coordinates.itertuples(index=False):
                factor_mean = represented_mean(
                    factor_bw, row.chrom, int(row.start), int(row.end)
                )
                igg_mean = represented_mean(
                    igg_bw, row.chrom, int(row.start), int(row.end)
                )
                ratio = (
                    factor_mean / igg_mean
                    if igg_mean > 0
                    else (math.inf if factor_mean > 0 else 0.0)
                )
                best = max(best, ratio)
            ratios.append(best)
        libraries = [
            f"{factor}-{index}"
            for index, ratio in enumerate(ratios, 1)
            if ratio >= RATIO_MIN
        ]
        rows.append(
            {
                "gene": gene,
                "promoter_factor_over_IgG_rep1": ratios[0],
                "promoter_factor_over_IgG_rep2": ratios[1],
                "n_libraries": len(libraries),
                "libraries": ";".join(libraries),
                "pooled_peak_pass": True,
            }
        )
    return pd.DataFrame(rows)


def write_bed(rows: list[tuple[str, int, int]], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(path, sep="\t", header=False, index=False)


def filtered_bam(library: str) -> Path:
    return CUTRUN_ROOT / "02_bam" / library / f"{library}.filt.sorted.bam"


def install_segments(
    segments: pd.DataFrame,
    promoter_sets: dict[str, set[str]],
    igg_peaks: list[tuple[str, int, int]],
) -> None:
    data = CUTRUN_ROOT / "data"
    loci = data / "PeakLoci"
    peak_sets = data / "PeakSets"
    figure = data / "figure_inputs"
    for directory in (
        loci,
        peak_sets,
        figure / "peak_loci",
        figure / "batch_stratified_peaks",
        figure / "peak_distribution",
    ):
        directory.mkdir(parents=True, exist_ok=True)
    canonical = segments[
        ["locus_id", "region_mask", "chrom", "start", "end", "membership"]
    ].copy()
    canonical.to_csv(loci / "Venn_Peaks_loci.tsv", sep="\t", index=False)
    canonical.to_csv(figure / "peak_loci" / "Venn_Peaks_loci.tsv", sep="\t", index=False)
    legacy = canonical[["chrom", "start", "end", "region_mask"]].copy()
    legacy["length_bp"] = legacy.end - legacy.start
    legacy["signature"] = legacy.region_mask.map(
        lambda mask: "".join(
            str(int(bool(int(mask) & bit))) for bit in (1, 2, 4)
        )
    )
    legacy["MCM3"] = legacy.region_mask.map(lambda mask: int(bool(int(mask) & 1)))
    legacy["NONO"] = legacy.region_mask.map(lambda mask: int(bool(int(mask) & 2)))
    legacy["PSPC1"] = legacy.region_mask.map(lambda mask: int(bool(int(mask) & 4)))
    legacy["segment_id"] = legacy.apply(
        lambda row: f"{row.chrom}:{int(row.start)}-{int(row.end)}", axis=1
    )
    legacy[
        [
            "chrom",
            "start",
            "end",
            "length_bp",
            "signature",
            "MCM3",
            "NONO",
            "PSPC1",
            "segment_id",
        ]
    ].to_csv(figure / "Venn_Peaks_segments.tsv", sep="\t", index=False)
    counts = (
        pd.DataFrame({"region_mask": range(1, 8)})
        .merge(
            canonical.groupby("region_mask").size().rename("peak_loci").reset_index(),
            how="left",
            on="region_mask",
        )
        .fillna({"peak_loci": 0})
    )
    counts["peak_loci"] = counts.peak_loci.astype(int)
    counts["minimum_peak_overlap_bp"] = PROMOTER_OVERLAP_BP
    counts["unit"] = "exact atomic genomic segment"
    for target in (data / "Venn_Peaks_counts.tsv", figure / "Venn_Peaks_counts.tsv"):
        counts.to_csv(target, sep="\t", index=False)
    summary = []
    for index, factor in enumerate(FACTORS):
        selected = canonical.loc[
            canonical.region_mask.map(
                lambda mask: bool(int(mask) & (1 << index))
            ),
            ["chrom", "start", "end"],
        ]
        for target in (
            loci / f"Venn_Peaks_{factor}.bed",
            peak_sets / f"Peaks_{factor}.bed",
            figure / "batch_stratified_peaks" / f"Venn_Peaks_{factor}.bed",
        ):
            selected.to_csv(target, sep="\t", header=False, index=False)
        summary.append(
            {
                "factor": factor,
                "peak_segments": len(selected),
                "promoter_bound_genes": len(promoter_sets[factor]),
                "rule": "pooled matched-IgG peak; promoter ratio >=2 in 2/2 replicates",
            }
        )
    write_bed(igg_peaks, peak_sets / "Peaks_IgG.bed")
    write_bed(igg_peaks, figure / "peak_distribution" / "Peaks_IgG.bed")
    pd.DataFrame(summary).to_csv(data / "Venn_Peaks_support.tsv", sep="\t", index=False)
    pd.DataFrame(summary).to_csv(figure / "Venn_Peaks_support.tsv", sep="\t", index=False)
    universe = [
        {"panel": panel, "set": factor, "unit": unit, "n": value}
        for factor in FACTORS
        for panel, unit, value in (
            (
                "Venn_Peaks",
                "factor-associated exact atomic genomic segments",
                int(summary[FACTORS.index(factor)]["peak_segments"]),
            ),
            (
                "Pie_PeakDistribution",
                "factor-associated exact atomic genomic segments",
                int(summary[FACTORS.index(factor)]["peak_segments"]),
            ),
            ("Venn_PromoterGenes", "promoter-bound genes", len(promoter_sets[factor])),
        )
    ]
    universe.append(
        {
            "panel": "Pie_IgG",
            "set": "IgG",
            "unit": "2/2 reproducible IgG intervals",
            "n": len(igg_peaks),
        }
    )
    pd.DataFrame(universe).to_csv(data / "Venn_Peaks_universe.tsv", sep="\t", index=False)


def install_promoters(
    promoter_sets: dict[str, set[str]],
    support_tables: dict[str, pd.DataFrame],
) -> None:
    data = CUTRUN_ROOT / "data"
    promoter_dir = CUTRUN_ROOT / "05_promoters"
    figure = data / "figure_inputs" / "promoter_gene_venn"
    promoter_dir.mkdir(parents=True, exist_ok=True)
    figure.mkdir(parents=True, exist_ok=True)
    for factor in FACTORS:
        genes = pd.DataFrame({"gene": sorted(promoter_sets[factor])})
        for target in (
            promoter_dir / f"Venn_PromoterGenes_{factor}.tsv",
            figure / f"Venn_PromoterGenes_{factor}.tsv",
        ):
            genes.to_csv(target, sep="\t", index=False)
        for target in (
            promoter_dir / f"Venn_PromoterGenes_{factor}_support.tsv",
            figure / f"Venn_PromoterGenes_{factor}_support.tsv",
        ):
            support_tables[factor].to_csv(target, sep="\t", index=False)
    regions = []
    for gene in sorted(set().union(*promoter_sets.values())):
        mask = sum(
            1 << index
            for index, factor in enumerate(FACTORS)
            if gene in promoter_sets[factor]
        )
        regions.append({"region": REGION_NAMES[mask], "gene": gene})
    pd.DataFrame(regions).to_csv(figure / "Venn_PromoterGenes_regions.tsv", sep="\t", index=False)
    summary = pd.DataFrame(
        [
            {
                "factor": factor,
                "promoter_bound_genes": len(promoter_sets[factor]),
                "rule": "pooled matched-IgG peak plus factor/IgG >=2 in 2/2 replicates",
            }
            for factor in FACTORS
        ]
    )
    summary.to_csv(data / "Venn_PromoterGenes_summary.tsv", sep="\t", index=False)
    parameters = pd.DataFrame(
        [
            ("analysis_role", "two-biological-replicate matched-IgG CUT&RUN analysis"),
            ("batches", "20200929 replicate 1; 20200923 replicate 2"),
            (
                "peak_call",
                "pooled two-replicate BAMPE against pooled matched IgG; candidate p<=1e-3",
            ),
            (
                "peak_thresholds",
                "retained p<=1e-4; q<=0.01; fold enrichment>=3",
            ),
            ("promoter_window", "GENCODE M25 TSS +/-1,000 bp"),
            ("peak_promoter_overlap", ">=250 bp"),
            (
                "binding_rule",
                "passing pooled peak plus promoter factor/IgG mean >=2 in 2/2 biological replicates",
            ),
            ("spike_in", "mouse paired-fragment coverage per 10,000 retained yeast pairs"),
        ],
        columns=["parameter", "value"],
    )
    parameters.to_csv(data / "AnalysisParameters.tsv", sep="\t", index=False)
    parameters.to_csv(figure / "AnalysisParameters.tsv", sep="\t", index=False)


def write_reference_beds(promoters: pd.DataFrame, bodies: pd.DataFrame) -> None:
    data = CUTRUN_ROOT / "data"
    data.mkdir(parents=True, exist_ok=True)
    promoters.to_csv(
        data / "Promoters_M25_TSSplusminus1kb.bed",
        sep="\t",
        header=False,
        index=False,
    )
    bodies.to_csv(data / "GeneBodies_M25.bed6", sep="\t", header=False, index=False)


def accepted() -> None:
    source = CUTRUN_ROOT / "data" / "figure_inputs"
    promoter_source = source / "promoter_gene_venn"
    promoter_sets, support_tables = {}, {}
    for factor in FACTORS:
        genes = pd.read_csv(
            promoter_source / f"Venn_PromoterGenes_{factor}.tsv",
            sep="\t",
            dtype=str,
        )
        promoter_sets[factor] = set(genes.gene.dropna())
        support_tables[factor] = pd.read_csv(
            promoter_source / f"Venn_PromoterGenes_{factor}_support.tsv",
            sep="\t",
        )
    raw_segments = pd.read_csv(source / "Venn_Peaks_segments.tsv", sep="\t")
    raw_segments["region_mask"] = (
        raw_segments.MCM3.astype(int)
        + 2 * raw_segments.NONO.astype(int)
        + 4 * raw_segments.PSPC1.astype(int)
    )
    raw_segments["locus_id"] = raw_segments.segment_id
    raw_segments["membership"] = raw_segments.region_mask.map(REGION_NAMES)
    igg_peaks = []
    with (
        CUTRUN_ROOT / "05_consensus" / "NT" / "IgG" / "union.merge.bed"
    ).open() as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 3:
                igg_peaks.append((fields[0], int(fields[1]), int(fields[2])))
    install_promoters(promoter_sets, support_tables)
    install_segments(raw_segments, promoter_sets, igg_peaks)
    promoters, bodies = gtf_tables(REFERENCE_DIR / "gencode.vM25.annotation.gtf")
    write_reference_beds(promoters, bodies)
    print("[DONE] Staged strict 2020 matched-IgG peak and promoter evidence.")


def raw(manifest: pd.DataFrame) -> None:
    macs3 = executable(
        "macs3", ("macs3", "/opt/anaconda3/envs/cutrun_env/bin/macs3")
    )
    promoters, bodies = gtf_tables(REFERENCE_DIR / "gencode.vM25.annotation.gtf")
    write_reference_beds(promoters, bodies)
    peak_root = CUTRUN_ROOT / "04_peaks"
    candidate_root = peak_root / "candidate_calls"
    peak_sets = {}
    promoter_sets = {}
    support_tables = {}
    igg_rows = manifest.loc[manifest.factor.eq("IgG")].sort_values("replicate")
    igg_bams = [
        filtered_bam(row.library_id) for row in igg_rows.itertuples(index=False)
    ]
    for factor in FACTORS:
        rows = manifest.loc[manifest.factor.eq(factor)].sort_values("replicate")
        bams = [filtered_bam(row.library_id) for row in rows.itertuples(index=False)]
        candidate = call_macs(
            macs3,
            bams,
            candidate_root / "pooled" / factor,
            f"{factor}-pooled_full",
            igg_bams,
        )
        peak = filter_peak_file(
            candidate,
            peak_root / "pooled" / factor / candidate.name,
        )
        peak_sets[factor] = merge_intervals(read_filtered_peaks(peak))
        candidates = overlap_genes(peak_sets[factor], promoters)
        support = promoter_support(promoters, candidates, factor)
        support_tables[factor] = support
        promoter_sets[factor] = set(support.loc[support.n_libraries.ge(2), "gene"])
    igg_calls = []
    for row in igg_rows.itertuples(index=False):
        candidate = call_macs(
            macs3,
            [filtered_bam(row.library_id)],
            candidate_root / f"rep{row.replicate}" / "IgG",
            f"{row.library_id}_full",
            None,
        )
        peak = filter_peak_file(
            candidate,
            peak_root / f"rep{row.replicate}" / "IgG" / candidate.name,
        )
        igg_calls.append(read_filtered_peaks(peak))
    igg_peaks = supported_segments(igg_calls, 2)
    install_promoters(promoter_sets, support_tables)
    install_segments(atomic_peak_segments(peak_sets), promoter_sets, igg_peaks)
    print("[DONE] Rebuilt strict 2020 matched-IgG peak and promoter evidence from BAMs.")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--from-raw", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()
    if args.dry_run:
        print("[DRY-RUN] Would define 2020 matched-IgG peaks and 2-of-2 promoter binding.")
        return
    manifest = pd.read_csv(CUTRUN_ROOT / "metadata" / "Samples.tsv", sep="\t")
    raw(manifest) if args.from_raw else accepted()


if __name__ == "__main__":
    main()
