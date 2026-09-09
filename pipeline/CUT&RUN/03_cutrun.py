#!/usr/bin/env python3
"""Generate final CUT&RUN alignments, filtered BAMs, and yeast-normalized tracks.

This stage processes the original 2019, 20200923, and 20200929 FASTQs.  It
does not consume BAMs or bigWigs from previous reanalysis directories.
"""

from __future__ import annotations

import argparse
import csv
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pandas as pd

PIPELINE_ROOT = Path(__file__).resolve().parents[1]
if str(PIPELINE_ROOT) not in sys.path:
    sys.path.insert(0, str(PIPELINE_ROOT))

from config import (
    CUTRUN_ROOT, DISPLAY_2020_TO_2019, DISPLAY_MULTIPLIER, FACTORS, MAPQ_MIN,
    REFERENCE_DIR, YEAST_SCALE_NUMERATOR,
)


def executable(name: str, candidates: tuple[str, ...]) -> str:
    for candidate in candidates:
        found = shutil.which(candidate)
        if found:
            return found
        if Path(candidate).is_file():
            return candidate
    raise FileNotFoundError(f"Required executable not found: {name}")


def run(command: list[str], dry_run: bool = False, stdout=None) -> None:
    print("[RUN]", " ".join(command), flush=True)
    if not dry_run:
        subprocess.run(command, check=True, stdout=stdout)


def align(bowtie2: str, samtools: str, index: Path, r1: Path, r2: Path, output: Path, threads: int, dry_run: bool) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    if output.is_file() and output.with_suffix(output.suffix + ".bai").is_file():
        return
    command1 = [bowtie2, "--local", "--very-sensitive-local", "--no-unal", "--no-mixed", "--no-discordant", "--phred33", "-I", "10", "-X", "700", "-x", str(index), "-1", str(r1), "-2", str(r2), "-p", str(threads)]
    command2 = [samtools, "sort", "-@", str(threads), "-o", str(output), "-"]
    print("[RUN]", " ".join(command1), "|", " ".join(command2), flush=True)
    if dry_run:
        return
    with subprocess.Popen(command1, stdout=subprocess.PIPE) as first:
        subprocess.run(command2, stdin=first.stdout, check=True)
        first.stdout.close()
        if first.wait() != 0:
            raise RuntimeError(f"Bowtie2 failed for {r1}")
    run([samtools, "index", "-@", str(threads), str(output)])


def mark_and_filter(picard: str, samtools: str, source: Path, destination: Path, threads: int, dry_run: bool) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.is_file() and destination.with_suffix(destination.suffix + ".bai").is_file():
        return
    marked = destination.with_name(destination.stem + ".markdup.bam")
    metrics = destination.with_name(destination.stem + ".duplication_metrics.txt")
    # Picard marks, rather than removes, duplicates.  The subsequent SAM flag
    # filter excludes duplicate-marked fragments identically for mouse/yeast.
    run([picard, "MarkDuplicates", f"I={source}", f"O={marked}", f"M={metrics}", "REMOVE_DUPLICATES=false", "ASSUME_SORTED=true", "VALIDATION_STRINGENCY=LENIENT"], dry_run)
    exclude = str(4 + 256 + 512 + 1024 + 2048)
    run([samtools, "view", "-@", str(threads), "-b", "-f", "2", "-F", exclude, "-q", str(MAPQ_MIN), "-o", str(destination), str(marked)], dry_run)
    run([samtools, "index", "-@", str(threads), str(destination)], dry_run)


def count_yeast_pairs(samtools: str, bam: Path) -> int:
    exclude = str(4 + 256 + 512 + 1024 + 2048)
    result = subprocess.run([samtools, "view", "-c", "-f", "66", "-F", exclude, "-q", str(MAPQ_MIN), str(bam)], check=True, capture_output=True, text=True)
    value = int(result.stdout.strip())
    if value <= 0:
        raise RuntimeError(f"No retained yeast read pairs in {bam}")
    return value


def make_track(bedtools: str, converter: str, bam: Path, chrom_sizes: Path, output: Path, scale: float, dry_run: bool) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    if output.is_file():
        return
    bedgraph = output.with_suffix(".bedGraph")
    command = [bedtools, "genomecov", "-bg", "-pc", "-ibam", str(bam), "-scale", f"{scale:.12g}"]
    print("[RUN]", " ".join(command), ">", bedgraph, flush=True)
    if dry_run:
        return
    with bedgraph.open("w") as handle:
        subprocess.run(command, check=True, stdout=handle)
    ordered = bedgraph.with_suffix(".sorted.bedGraph")
    with ordered.open("w") as handle:
        subprocess.run(["sort", "-k1,1", "-k2,2n", str(bedgraph)], check=True, stdout=handle)
    run([converter, str(ordered), str(chrom_sizes), str(output)])
    bedgraph.unlink(missing_ok=True); ordered.unlink(missing_ok=True)


def mean_and_scale_tracks(manifest: pd.DataFrame) -> None:
    """Materialize exact equal-weight batch means from final per-library tracks."""
    try:
        import pyBigWig
    except ImportError as error:
        raise RuntimeError("pyBigWig is required for replicate mean tracks") from error

    track_root = CUTRUN_ROOT / "03_bigwig"
    means = track_root / "mean_intermediates"
    igv = track_root / "IGV_representation"
    means.mkdir(parents=True, exist_ok=True); igv.mkdir(parents=True, exist_ok=True)

    def scale_bigwig(source: Path, destination: Path, multiplier: float) -> None:
        if destination.is_file():
            return
        with pyBigWig.open(str(source)) as reader, pyBigWig.open(str(destination), "w") as writer:
            chroms = reader.chroms(); writer.addHeader(list(chroms.items()), maxZooms=10)
            for chrom in chroms:
                intervals = reader.intervals(chrom) or []
                if intervals:
                    writer.addEntries([chrom] * len(intervals), [entry[0] for entry in intervals], ends=[entry[1] for entry in intervals], values=[entry[2] * multiplier for entry in intervals])

    # The final track means retain every original coverage interval.  The two
    # batches are averaged with equal batch weight, not unequal 4:2 library weight.
    def mean_bigwig(sources: list[Path], destination: Path) -> None:
        if destination.is_file():
            return
        readers = [pyBigWig.open(str(path)) for path in sources]
        try:
            chroms = readers[0].chroms()
            with pyBigWig.open(str(destination), "w") as writer:
                writer.addHeader(list(chroms.items()), maxZooms=10)
                for chrom, length in chroms.items():
                    boundaries = {0, int(length)}
                    interval_sets = []
                    for reader in readers:
                        rows = reader.intervals(chrom) or []
                        interval_sets.append(rows)
                        for start, end, _ in rows:
                            boundaries.add(int(start)); boundaries.add(int(end))
                    points = sorted(boundaries)
                    starts: list[int] = []; ends: list[int] = []; values: list[float] = []
                    positions = [0] * len(readers)
                    for left, right in zip(points[:-1], points[1:]):
                        value = 0.0
                        for idx, rows in enumerate(interval_sets):
                            while positions[idx] < len(rows) and rows[positions[idx]][1] <= left:
                                positions[idx] += 1
                            if positions[idx] < len(rows):
                                start, end, coverage = rows[positions[idx]]
                                if start <= left and end >= right:
                                    value += float(coverage)
                        value /= len(readers)
                        if value != 0.0:
                            starts.append(left); ends.append(right); values.append(value)
                    if starts:
                        writer.addEntries([chrom] * len(starts), starts, ends=ends, values=values)
        finally:
            for reader in readers:
                reader.close()

    rows = []
    for factor in FACTORS:
        tracks19 = [track_root / str(row.library_id) / f"{row.library_id}.bw" for row in manifest.itertuples() if row.batch == "2019" and bool(row.selected_factor_library) and row.factor == factor]
        tracks20 = [track_root / str(row.library_id) / f"{row.library_id}.bw" for row in manifest.itertuples() if row.batch == "2020" and row.factor == factor]
        if len(tracks19) != 4 or len(tracks20) != 2:
            raise RuntimeError(f"Expected 4 2019 and 2 2020 tracks for {factor}; got {len(tracks19)} and {len(tracks20)}")
        mean19 = means / f"{factor}_2019_mean.bw"; mean20 = means / f"{factor}_2020_mean.bw"
        calibrated20 = means / f"{factor}_2020_to_2019.bw"; combined = means / f"{factor}_six_library_batchBalanced_mean.bw"
        mean_bigwig(tracks19, mean19); mean_bigwig(tracks20, mean20); scale_bigwig(mean20, calibrated20, DISPLAY_2020_TO_2019[factor]); mean_bigwig([mean19, calibrated20], combined)
        scale_bigwig(combined, igv / f"{factor}_mean.bw", DISPLAY_MULTIPLIER)
        rows.append({"factor": factor, "n_2019": 4, "n_2020": 2, "2020_to_2019_display_multiplier": DISPLAY_2020_TO_2019[factor]})
    igg = [track_root / str(row.library_id) / f"{row.library_id}.bw" for row in manifest.itertuples() if row.batch == "2020" and row.factor == "IgG"]
    if len(igg) != 2:
        raise RuntimeError("Expected two 2020 IgG tracks")
    igg_mean = means / "IgG_2020_mean.bw"; igg_calibrated = means / "IgG_2020_to_2019_commonBridge_mean.bw"
    mean_bigwig(igg, igg_mean); scale_bigwig(igg_mean, igg_calibrated, DISPLAY_2020_TO_2019["IgG"]); scale_bigwig(igg_calibrated, igv / "IgG_mean.bw", DISPLAY_MULTIPLIER)
    pd.DataFrame(rows).to_csv(means / "BatchBalancedTrackSummary.tsv", sep="\t", index=False)
    source_tracks = [f"{factor}_six_library_batchBalanced_mean.bw" for factor in FACTORS]
    source_tracks.append("IgG_2020_to_2019_commonBridge_mean.bw")
    igv_tracks = [f"{factor}_mean.bw" for factor in FACTORS]
    igv_tracks.append("IgG_mean.bw")
    pd.DataFrame({
        "source_track": [f"mean_intermediates/{name}" for name in source_tracks],
        "igv_track": igv_tracks,
        "display_multiplier": DISPLAY_MULTIPLIER,
        "analysis_role": "publication visualization only",
    }).to_csv(igv / "IGV_tracks_manifest.tsv", sep="\t", index=False)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--threads", type=int, default=max(1, min(16, os.cpu_count() or 1)))
    parser.add_argument("--from-raw", action="store_true", help="Explicitly rebuild tracks from original FASTQs")
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()
    if not args.from_raw:
        raise SystemExit(
            "03_cutrun.py is the raw regeneration path. Use run_all.py --source raw, "
            "or the default accepted-intermediate workflow for figure production."
        )
    manifest_path = CUTRUN_ROOT / "metadata" / "Samples.tsv"
    if not manifest_path.is_file():
        raise FileNotFoundError("Run 01_prepare.py before CUT&RUN processing")
    manifest = pd.read_csv(manifest_path, sep="\t")
    bowtie2 = executable("bowtie2", ("bowtie2", "/opt/anaconda3/envs/cutrun_env/bin/bowtie2"))
    samtools = executable("samtools", ("samtools", "/opt/anaconda3/envs/cutrun_env/bin/samtools"))
    picard = executable("picard", ("picard", "/opt/anaconda3/envs/cutrun_env/bin/picard"))
    bedtools = executable("bedtools", ("bedtools", "/opt/anaconda3/envs/cutrun_env/bin/bedtools"))
    converter = executable("bedGraphToBigWig", ("bedGraphToBigWig", "/opt/anaconda3/envs/cutrun_env/bin/bedGraphToBigWig"))
    mouse_index = REFERENCE_DIR / "bowtie2_mm10" / "mm10"
    yeast_index = REFERENCE_DIR / "bowtie2_sacCer3" / "sacCer3"
    chrom_sizes = REFERENCE_DIR / "mm10.chrom.sizes"
    if not args.dry_run:
        index_requirements = ((mouse_index, "mm10 Bowtie2 index"), (yeast_index, "sacCer3 Bowtie2 index"))
        for prefix, label in index_requirements:
            if not list(prefix.parent.glob(prefix.name + ".*.bt2*")):
                raise FileNotFoundError(f"Missing {label}: {prefix}.{{1..4}}.bt2")
        if not chrom_sizes.is_file():
            raise FileNotFoundError(f"Missing mm10 chromosome sizes: {chrom_sizes}")

    selected = manifest.loc[
        (manifest["batch"].astype(str) == "2020")
        | manifest["selected_factor_library"].astype(str).str.lower().eq("true")
    ].copy()
    rows: list[dict[str, object]] = []
    for row in selected.itertuples(index=False):
        library = str(row.library_id); r1 = Path(str(row.r1)); r2 = Path(str(row.r2))
        mouse_raw = CUTRUN_ROOT / "01_bowtie2" / "mouse" / f"{library}.bam"
        yeast_raw = CUTRUN_ROOT / "01_bowtie2" / "yeast" / f"{library}.bam"
        mouse_filtered = CUTRUN_ROOT / "02_bam" / f"{library}.bam"
        yeast_filtered = CUTRUN_ROOT / "02_bam" / "yeast" / f"{library}.bam"
        align(bowtie2, samtools, mouse_index, r1, r2, mouse_raw, args.threads, args.dry_run)
        align(bowtie2, samtools, yeast_index, r1, r2, yeast_raw, args.threads, args.dry_run)
        mark_and_filter(picard, samtools, mouse_raw, mouse_filtered, args.threads, args.dry_run)
        mark_and_filter(picard, samtools, yeast_raw, yeast_filtered, args.threads, args.dry_run)
        if args.dry_run:
            continue
        yeast_pairs = count_yeast_pairs(samtools, yeast_filtered)
        scale = YEAST_SCALE_NUMERATOR / yeast_pairs
        bw = CUTRUN_ROOT / "03_bigwig" / library / f"{library}.bw"
        make_track(bedtools, converter, mouse_filtered, chrom_sizes, bw, scale, args.dry_run)
        rows.append({"library_id": library, "yeast_pairs": yeast_pairs, "coverage_scale": scale, "normalization": "mouse paired-fragment coverage per 10,000 retained yeast pairs"})
    if not args.dry_run:
        pd.DataFrame(rows).to_csv(CUTRUN_ROOT / "03_bigwig" / "YeastNormalization.tsv", sep="\t", index=False)
        mean_and_scale_tracks(manifest)


if __name__ == "__main__":
    main()
