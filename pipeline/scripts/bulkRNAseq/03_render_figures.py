#!/usr/bin/env python3
"""Render all bulk RNA-seq publication figures from packaged quantifications."""

from __future__ import annotations

import argparse
import gzip
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pandas as pd


SCRIPT_DIR = Path(__file__).resolve().parent
PIPELINE_ROOT = SCRIPT_DIR.parent
if str(PIPELINE_ROOT) not in sys.path:
    sys.path.insert(0, str(PIPELINE_ROOT))

from config import FACTORS, MM10_FASTA_NAME, REFERENCE_DIR, RNA_FDR_MAX, RNA_LOGFC_MIN, RNA_ROOT


VENN_RENDERER = SCRIPT_DIR / "figures_rendering" / "venn.py"
OUTPUT = RNA_ROOT / "visuals"
PUBLICATION_FIGURES = OUTPUT / "publication_figures"
TRACK_GROUPS = {
    "NT": lambda row: row["condition"] == "NT",
    "shMCM3": lambda row: row["condition"] == "KD" and row["target"] == "MCM3",
    "shNONO": lambda row: row["condition"] == "KD" and row["target"] == "NONO",
    "shPSPC1": lambda row: row["condition"] == "KD" and row["target"] == "PSPC1",
}


def executable(name: str) -> str:
    path = shutil.which(name)
    if path:
        return path
    environment_path = Path(sys.prefix) / "bin" / name
    if environment_path.is_file() and environment_path.stat().st_mode & 0o111:
        return str(environment_path)
    raise FileNotFoundError(f"Required executable not found on PATH: {name}")


def run(command: list[str], dry_run: bool = False) -> None:
    print("[RUN]", " ".join(command), flush=True)
    if not dry_run:
        subprocess.run(command, check=True)


def selected_samples() -> pd.DataFrame:
    samples = pd.read_csv(RNA_ROOT / "metadata" / "Samples.tsv", sep="\t")
    samples = samples.loc[samples["selected"].astype(str).str.lower().eq("true")].copy()
    if samples.empty:
        raise ValueError("No selected RNA-seq samples found in metadata/Samples.tsv")
    return samples


def read_length(fastq: Path) -> int:
    with gzip.open(fastq, "rt") as handle:
        next(handle)
        return len(next(handle).strip())


def star_index(star: str, samples: pd.DataFrame, output: Path, threads: int, dry_run: bool) -> Path:
    index = output / "reference" / "STAR_GRCm38_GENCODE_M25"
    if (index / "Genome").is_file():
        return index
    fasta = REFERENCE_DIR / MM10_FASTA_NAME
    gtf = REFERENCE_DIR / "gencode.vM25.annotation.gtf"
    if not fasta.is_file() or not gtf.is_file():
        raise FileNotFoundError("Missing mm10 FASTA or GENCODE M25 GTF in pipeline/reference")
    index.mkdir(parents=True, exist_ok=True)
    run(
        [
            star, "--runMode", "genomeGenerate", "--runThreadN", str(threads),
            "--genomeDir", str(index), "--genomeFastaFiles", str(fasta),
            "--sjdbGTFfile", str(gtf), "--sjdbOverhang",
            str(max(read_length(Path(path)) for path in samples["r1"]) - 1),
            "--genomeSAindexNbases", "13",
        ],
        dry_run,
    )
    return index


def align_sample(star: str, samtools: str, bam_coverage: str, index: Path, row: pd.Series, output: Path, threads: int, dry_run: bool) -> tuple[Path, Path]:
    sample = str(row["sample_id"])
    sample_dir = output / "alignments" / sample
    sample_dir.mkdir(parents=True, exist_ok=True)
    temporary_dir = output / ".tmp" / sample
    sorted_bam = sample_dir / "Aligned.sortedByCoord.out.bam"
    filtered_bam = sample_dir / f"{sample}.primary.proper.MAPQ30.bam"
    track = output / "individual" / f"{sample}.bw"
    track.parent.mkdir(parents=True, exist_ok=True)
    if not filtered_bam.is_file():
        if not sorted_bam.is_file():
            if not dry_run and temporary_dir.exists():
                shutil.rmtree(temporary_dir)
            try:
                run(
                    [
                        star, "--runThreadN", str(threads), "--genomeDir", str(index),
                        "--readFilesIn", str(row["r1"]), str(row["r2"]), "--readFilesCommand", "zcat",
                        "--twopassMode", "Basic", "--outSAMtype", "BAM", "SortedByCoordinate",
                        "--outSAMattributes", "NH", "HI", "AS", "nM", "XS",
                        "--outSAMattrRGline", f"ID:{sample}", f"SM:{sample}", "PL:ILLUMINA",
                        "--outFileNamePrefix", f"{sample_dir}/", "--outTmpDir", str(temporary_dir),
                    ],
                    dry_run,
                )
            finally:
                if not dry_run and temporary_dir.exists():
                    shutil.rmtree(temporary_dir)
        run([samtools, "view", "-@", str(threads), "-b", "-q", "30", "-f", "2", "-F", "2304", "-o", str(filtered_bam), str(sorted_bam)], dry_run)
        run([samtools, "index", "-@", str(threads), str(filtered_bam)], dry_run)
        if not dry_run and sorted_bam.is_file():
            sorted_bam.unlink()
    if not track.is_file():
        run([bam_coverage, "--bam", str(filtered_bam), "--outFileName", str(track), "--outFileFormat", "bigwig", "--normalizeUsing", "CPM", "--binSize", "10", "--numberOfProcessors", str(threads)], dry_run)
    return filtered_bam, track


def mean_track(bigwig_compare: str, tracks: list[Path], output: Path, threads: int, dry_run: bool) -> None:
    if output.is_file():
        return
    if not tracks:
        raise ValueError(f"Cannot average zero tracks for {output.name}")
    output.parent.mkdir(parents=True, exist_ok=True)
    working = output.with_suffix(".working.bw")
    if dry_run:
        print(f"[RUN] copy {tracks[0]} {working}")
    else:
        shutil.copy2(tracks[0], working)
    for number, track in enumerate(tracks[1:], start=2):
        next_working = output.with_suffix(f".working_{number}.bw")
        run([bigwig_compare, "--bigwig1", str(working), "--bigwig2", str(track), "--operation", "add", "--scaleFactors", f"{(number - 1) / number}:{1 / number}", "--binSize", "10", "--numberOfProcessors", str(threads), "--outFileName", str(next_working)], dry_run)
        if not dry_run:
            working.unlink()
            next_working.replace(working)
    if not dry_run:
        working.replace(output)


def generate_igv_tracks(threads: int, dry_run: bool) -> None:
    output = RNA_ROOT / "igv_tracks"
    samples = selected_samples()
    index = star_index(executable("STAR"), samples, output, threads, dry_run)
    samtools = executable("samtools")
    bam_coverage = executable("bamCoverage")
    tracks: dict[str, Path] = {}
    records: list[dict[str, str]] = []
    for _, row in samples.iterrows():
        bam, track = align_sample(executable("STAR"), samtools, bam_coverage, index, row, output, threads, dry_run)
        sample = str(row["sample_id"])
        tracks[sample] = track
        records.append({"sample_id": sample, "group": next(group for group, predicate in TRACK_GROUPS.items() if predicate(row)), "experiment": str(row["experiment"]), "condition": str(row["condition"]), "target": str(row["target"]), "shRNA": str(row["shRNA"]), "bam": str(bam), "bigwig": str(track), "normalization": "CPM; 10-bp bins; primary properly paired MAPQ>=30 alignments"})
    for group, predicate in TRACK_GROUPS.items():
        ids = samples.loc[samples.apply(predicate, axis=1), "sample_id"].astype(str).tolist()
        mean_track(executable("bigwigCompare"), [tracks[sample] for sample in ids], output / "means" / f"{group}_mean.bw", threads, dry_run)
    if not dry_run:
        pd.DataFrame(records).to_csv(output / "IGV_tracks_manifest.tsv", sep="\t", index=False)
    print(f"[DONE] RNA-seq IGV tracks: {output / 'means'}")


def genes_for_direction(factor: str, direction: str) -> set[str]:
    table = pd.read_csv(
        RNA_ROOT / "data" / f"Volcano_{factor}_data.tsv",
        sep="\t",
        dtype={"gene_id": str},
    )
    keep = (pd.to_numeric(table["adj.P.Val"], errors="coerce") <= RNA_FDR_MAX) & (
        pd.to_numeric(table["logFC"], errors="coerce").abs() >= RNA_LOGFC_MIN
    )
    logfc = pd.to_numeric(table["logFC"], errors="coerce")
    keep &= logfc.gt(0) if direction == "UP" else logfc.lt(0)
    return set(table.loc[keep, "gene_id"].dropna().astype(str))


def venn_regions(sets: dict[str, set[str]]) -> dict[str, int]:
    mcm3, nono, pspc1 = (sets[factor] for factor in FACTORS)
    return {
        "n100": len(mcm3 - nono - pspc1),
        "n010": len(nono - mcm3 - pspc1),
        "n001": len(pspc1 - mcm3 - nono),
        "n110": len((mcm3 & nono) - pspc1),
        "n101": len((mcm3 & pspc1) - nono),
        "n011": len((nono & pspc1) - mcm3),
        "n111": len(mcm3 & nono & pspc1),
    }


def render_venn(direction: str, output: Path, no_text: bool = False) -> None:
    if not VENN_RENDERER.is_file():
        raise FileNotFoundError(VENN_RENDERER)
    sets = {factor: genes_for_direction(factor, direction) for factor in FACTORS}
    title = "Unregulated genes" if direction == "UP" else "Downregulated genes"
    output.mkdir(parents=True, exist_ok=True)
    environment = os.environ.copy()
    environment.setdefault("MPLCONFIGDIR", str(RNA_ROOT.parent / ".cache" / "matplotlib"))
    counts = venn_regions(sets)
    for suffix, hide_numbers in (("_noTexts", True),) if no_text else (("", False),):
        command = [
            sys.executable,
            str(VENN_RENDERER),
            "--out",
            str(output / f"VennDiagram_{direction}{suffix}.png"),
            "--title",
            title,
            "--a-name",
            "MCM3",
            "--b-name",
            "NONO",
            "--c-name",
            "PSPC1",
            "--a-total",
            str(len(sets["MCM3"])),
            "--b-total",
            str(len(sets["NONO"])),
            "--c-total",
            str(len(sets["PSPC1"])),
            *[item for pair in counts.items() for item in (f"--{pair[0]}", str(pair[1]))],
        ]
        if hide_numbers:
            command.append("--hide-numbers")
        subprocess.run(command, check=True, env=environment)
    pd.DataFrame(
        [
            {"factor": factor, "direction": direction, "gene_id": gene}
            for factor, genes in sets.items()
            for gene in sorted(genes)
        ]
    ).to_csv(RNA_ROOT / "data" / f"VennDiagram_{direction}_data.tsv", sep="\t", index=False)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--publication-figures", action="store_true", help="Also render text-free PNGs in deg_work/visuals/publication_figures")
    parser.add_argument("--igv-tracks", action="store_true", help="Also generate CPM-normalized RNA-seq IGV tracks from selected raw libraries")
    args = parser.parse_args()
    if args.threads < 1:
        raise ValueError("--threads must be positive")

    manifest = RNA_ROOT / "metadata" / "Salmon_quantifications.tsv"
    if not manifest.is_file() and not args.dry_run:
        raise FileNotFoundError(f"Missing packaged quantification manifest: {manifest}")

    r_command = [
        executable("Rscript"),
        str(SCRIPT_DIR / "figures_rendering" / "render_deg_figures.R"),
        str(RNA_ROOT),
        str(REFERENCE_DIR),
    ]
    if args.publication_figures:
        r_command.append(str(PUBLICATION_FIGURES))
    if args.dry_run:
        print("[DRY-RUN]", " ".join(r_command))
        print("[DRY-RUN] Would render RNA-seq UP and DOWN Venn diagrams.")
        if args.igv_tracks:
            generate_igv_tracks(args.threads, dry_run=True)
        return

    subprocess.run(r_command, check=True)
    for direction in ("UP", "DOWN"):
        render_venn(direction, OUTPUT)
    if args.publication_figures:
        for direction in ("UP", "DOWN"):
            render_venn(direction, PUBLICATION_FIGURES, no_text=True)
        print(f"[DONE] Bulk RNA-seq publication figures: {PUBLICATION_FIGURES}")
    print(f"[DONE] Bulk RNA-seq figures: {RNA_ROOT / 'visuals'}")
    if args.igv_tracks:
        generate_igv_tracks(args.threads, dry_run=False)


if __name__ == "__main__":
    main()
