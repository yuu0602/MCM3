#!/usr/bin/env python3
"""Align the CUT&RUN experiments and build normalized tracks."""

#Before you run this script, please replace directory placeholders with your own directory

from __future__ import annotations

import argparse
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
    CUTRUN_ROOT,
    DISPLAY_MULTIPLIER,
    FACTORS,
    MAPQ_MIN,
    REFERENCE_DIR,
    YEAST_SCALE_NUMERATOR,
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


def indexed_bam_exists(path: Path) -> bool:
    return path.is_file() and path.with_suffix(path.suffix + ".bai").is_file()


def align(
    bowtie2: str,
    samtools: str,
    index: Path,
    r1: Path,
    r2: Path,
    output: Path,
    threads: int,
    yeast: bool,
    dry_run: bool,
) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    if indexed_bam_exists(output):
        return
    sensitivity = ["--very-sensitive", "-k", "2"] if yeast else ["--local", "--very-sensitive-local"]
    first = [
        bowtie2,
        *sensitivity,
        "--no-unal",
        "--no-mixed",
        "--no-discordant",
        "--phred33",
        "-I",
        "10",
        "-X",
        "700",
        "-x",
        str(index),
        "-1",
        str(r1),
        "-2",
        str(r2),
        "-p",
        str(threads),
    ]
    second = [samtools, "sort", "-@", str(threads), "-o", str(output), "-"]
    print("[RUN]", " ".join(first), "|", " ".join(second), flush=True)
    if dry_run:
        return
    with subprocess.Popen(first, stdout=subprocess.PIPE) as process:
        subprocess.run(second, stdin=process.stdout, check=True)
        process.stdout.close()
        if process.wait() != 0:
            raise RuntimeError(f"Bowtie2 failed for {r1}")
    run([samtools, "index", "-@", str(threads), str(output)])


def mark_and_filter(
    picard: str,
    samtools: str,
    source: Path,
    destination: Path,
    library: str,
    threads: int,
    dry_run: bool,
    include_duplicates: bool,
) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    if indexed_bam_exists(destination):
        return
    if include_duplicates:
        # Marking only changes the duplicate flag. This branch retains those reads,
        # so filter the original coordinate-sorted alignment directly.
        run(
            [
                samtools,
                "view",
                "-@",
                str(threads),
                "-b",
                "-f",
                "3",
                "-F",
                "2816",
                "-q",
                str(MAPQ_MIN),
                "-o",
                str(destination),
                str(source),
            ],
            dry_run,
        )
        run([samtools, "index", "-@", str(threads), str(destination)], dry_run)
        return
    read_groups = destination.with_name(f"{library}.rg.bam")
    marked = destination.with_name(f"{library}.markdup.bam")
    metrics = destination.with_name("markdup.metrics.txt")
    # 3840 excludes secondary, QC-failed, duplicate-marked, and supplementary reads.
    excluded_flags = "3840"
    run(
        [
            picard,
            "AddOrReplaceReadGroups",
            f"I={source}",
            f"O={read_groups}",
            "RGID=1",
            f"RGLB={library}",
            "RGPL=ILLUMINA",
            f"RGPU={library}",
            f"RGSM={library}",
            "VALIDATION_STRINGENCY=LENIENT",
        ],
        dry_run,
    )
    run(
        [
            picard,
            "MarkDuplicates",
            f"I={read_groups}",
            f"O={marked}",
            f"M={metrics}",
            "REMOVE_DUPLICATES=false",
            "ASSUME_SORTED=true",
            "VALIDATION_STRINGENCY=LENIENT",
        ],
        dry_run,
    )
    run(
        [
            samtools,
            "view",
            "-@",
            str(threads),
            "-b",
            "-f",
            "3",
            "-F",
            excluded_flags,
            "-q",
            str(MAPQ_MIN),
            "-o",
            str(destination),
            str(marked),
        ],
        dry_run,
    )
    run([samtools, "index", "-@", str(threads), str(destination)], dry_run)


def count_yeast_pairs(samtools: str, bam: Path) -> int:
    result = subprocess.run(
        [
            samtools,
            "view",
            "-c",
            "-f",
            "67",
            "-F",
            "2304",
            "-q",
            str(MAPQ_MIN),
            str(bam),
        ],
        check=True,
        capture_output=True,
        text=True,
    )
    count = int(result.stdout.strip())
    if count <= 0:
        raise RuntimeError(f"No retained yeast read pairs in {bam}")
    return count


def make_track(
    bedtools: str,
    converter: str,
    bam: Path,
    chrom_sizes: Path,
    output: Path,
    scale: float,
    dry_run: bool,
) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    if output.is_file():
        return
    bedgraph = output.with_suffix(".bedGraph")
    ordered = output.with_suffix(".sorted.bedGraph")
    command = [
        bedtools,
        "genomecov",
        "-bg",
        "-pc",
        "-ibam",
        str(bam),
        "-scale",
        f"{scale:.12g}",
    ]
    print("[RUN]", " ".join(command), ">", bedgraph, flush=True)
    if dry_run:
        return
    with bedgraph.open("w") as handle:
        subprocess.run(command, check=True, stdout=handle)
    with ordered.open("w") as handle:
        subprocess.run(["sort", "-k1,1", "-k2,2n", str(bedgraph)], check=True, stdout=handle)
    run([converter, str(ordered), str(chrom_sizes), str(output)])
    bedgraph.unlink(missing_ok=True)
    ordered.unlink(missing_ok=True)


def transform_bigwigs(sources: list[Path], destination: Path, multiplier: float = 1.0) -> None:
    import pyBigWig

    readers = [pyBigWig.open(str(path)) for path in sources]
    try:
        reader_chroms = [reader.chroms() for reader in readers]
        chroms: dict[str, int] = {}
        for header in reader_chroms:
            for chrom, length in header.items():
                known = chroms.setdefault(chrom, int(length))
                if known != int(length):
                    raise RuntimeError(f"Inconsistent chromosome length for {chrom}")
        destination.parent.mkdir(parents=True, exist_ok=True)
        with pyBigWig.open(str(destination), "w") as writer:
            writer.addHeader(list(chroms.items()), maxZooms=10)
            for chrom, length in chroms.items():
                boundaries = {0, int(length)}
                rows_by_reader = []
                for reader, header in zip(readers, reader_chroms):
                    rows = reader.intervals(chrom, 0, int(length)) if chrom in header else []
                    rows_by_reader.append(rows)
                    for start, end, _ in rows:
                        boundaries.update((int(start), int(end)))
                positions = [0] * len(readers)
                starts = []
                ends = []
                values = []
                points = sorted(boundaries)
                for left, right in zip(points[:-1], points[1:]):
                    total = 0.0
                    for index, rows in enumerate(rows_by_reader):
                        while positions[index] < len(rows) and rows[positions[index]][1] <= left:
                            positions[index] += 1
                        if positions[index] < len(rows):
                            start, end, value = rows[positions[index]]
                            if start <= left and end >= right:
                                total += float(value)
                    value = total / len(readers) * multiplier
                    if value:
                        starts.append(left)
                        ends.append(right)
                        values.append(value)
                if starts:
                    writer.addEntries([chrom] * len(starts), starts, ends=ends, values=values)
    finally:
        for reader in readers:
            reader.close()


def make_mean_tracks(manifest: pd.DataFrame, cutrun_root: Path) -> None:
    root = cutrun_root / "03_bigwig"
    base = root / "mean_intermediates"
    display = root / "IGV_representation"
    rows = []
    for factor in (*FACTORS, "IgG"):
        selected = manifest.loc[manifest.factor.eq(factor)].sort_values("replicate")
        if len(selected) != 2:
            raise RuntimeError(f"Expected two {factor} libraries; found {len(selected)}")
        tracks = [
            root / row.library_id / f"{row.library_id}.bw"
            for row in selected.itertuples(index=False)
        ]
        mean = base / f"{factor}_mean.bw"
        transform_bigwigs(tracks, mean)
        transform_bigwigs([mean], display / f"{factor}_mean.bw", DISPLAY_MULTIPLIER)
        rows.append(
            {
                "factor": factor,
                "replicates": 2,
                "source_track": str(mean.relative_to(root)),
                "igv_track": f"{factor}_mean.bw",
                "display_multiplier": DISPLAY_MULTIPLIER,
            }
        )
    pd.DataFrame(rows).to_csv(display / "IGV_tracks_manifest.tsv", sep="\t", index=False)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--threads", type=int, default=max(1, min(16, os.cpu_count() or 1))
    )
    parser.add_argument("--from-raw", action="store_true")
    parser.add_argument(
        "--work-root",
        type=Path,
        default=CUTRUN_ROOT,
        help="CUT&RUN output root; defaults to the primary analysis root.",
    )
    parser.add_argument(
        "--reuse-alignments",
        type=Path,
        help="Reuse mouse and yeast coordinate-sorted BAMs from this 01_bowtie2 directory.",
    )
    parser.add_argument(
        "--exclude-duplicates",
        dest="include_duplicates",
        action="store_false",
        help="Exclude duplicate-marked reads after Picard marking.",
    )
    parser.set_defaults(include_duplicates=True)
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()
    if not args.from_raw:
        raise SystemExit("Use run_pipeline.py --source raw to rebuild CUT&RUN intermediates")
    cutrun_root = args.work_root.resolve()
    manifest = pd.read_csv(CUTRUN_ROOT / "metadata" / "Samples.tsv", sep="\t")
    if len(manifest) != 8 or set(manifest.factor) != {*FACTORS, "IgG"}:
        raise RuntimeError("The CUT&RUN manifest must contain two libraries for each factor and IgG")
    bowtie2 = executable(
        "bowtie2", ("bowtie2", "/opt/anaconda3/envs/cutrun_env/bin/bowtie2")
    )
    samtools = executable(
        "samtools", ("samtools", "/opt/anaconda3/envs/cutrun_env/bin/samtools")
    )
    picard = executable(
        "picard", ("picard", "/opt/anaconda3/envs/cutrun_env/bin/picard")
    )
    bedtools = executable(
        "bedtools", ("bedtools", "/opt/anaconda3/envs/cutrun_env/bin/bedtools")
    )
    converter = executable(
        "bedGraphToBigWig",
        ("bedGraphToBigWig", "/opt/anaconda3/envs/cutrun_env/bin/bedGraphToBigWig"),
    )
    mouse_index = REFERENCE_DIR / "bowtie2_mm10" / "mm10"
    yeast_index = REFERENCE_DIR / "bowtie2_sacCer3" / "sacCer3"
    chrom_sizes = REFERENCE_DIR / "mm10.chrom.sizes"
    rows = []
    for sample in manifest.itertuples(index=False):
        library = str(sample.library_id)
        alignment_root = args.reuse_alignments or (cutrun_root / "01_bowtie2")
        mouse_raw = alignment_root / library / f"{library}.sorted.bam"
        yeast_raw = alignment_root / "yeast" / library / f"{library}.sorted.bam"
        mouse = cutrun_root / "02_bam" / library / f"{library}.filt.sorted.bam"
        yeast = cutrun_root / "02_bam" / "yeast" / library / f"{library}.filt.sorted.bam"
        if args.reuse_alignments:
            if not indexed_bam_exists(mouse_raw) or not indexed_bam_exists(yeast_raw):
                raise FileNotFoundError(f"Missing reusable alignment for {library}")
        else:
            align(
                bowtie2,
                samtools,
                mouse_index,
                Path(sample.r1),
                Path(sample.r2),
                mouse_raw,
                args.threads,
                False,
                args.dry_run,
            )
            align(
                bowtie2,
                samtools,
                yeast_index,
                Path(sample.r1),
                Path(sample.r2),
                yeast_raw,
                args.threads,
                True,
                args.dry_run,
            )
        mark_and_filter(
            picard, samtools, mouse_raw, mouse, library, args.threads, args.dry_run,
            args.include_duplicates,
        )
        mark_and_filter(
            picard, samtools, yeast_raw, yeast, library, args.threads, args.dry_run,
            args.include_duplicates,
        )
        if args.dry_run:
            continue
        yeast_pairs = count_yeast_pairs(samtools, yeast)
        scale = YEAST_SCALE_NUMERATOR / yeast_pairs
        make_track(
            bedtools,
            converter,
            mouse,
            chrom_sizes,
            cutrun_root / "03_bigwig" / library / f"{library}.bw",
            scale,
            False,
        )
        rows.append(
            {
                "library_id": library,
                "yeast_pairs": yeast_pairs,
                "coverage_scale": scale,
                "normalization": "mouse paired-fragment coverage per 10,000 retained yeast pairs",
                "duplicate_policy": "included" if args.include_duplicates else "excluded",
            }
        )
    if not args.dry_run:
        pd.DataFrame(rows).to_csv(cutrun_root / "03_bigwig" / "YeastNormalization.tsv", sep="\t", index=False)
        make_mean_tracks(manifest, cutrun_root)


if __name__ == "__main__":
    main()
