#!/usr/bin/env python3
"""Download/pin GENCODE M25; optionally build Salmon/Bowtie2 indexes."""

from __future__ import annotations

import argparse
import gzip
import shutil
import subprocess
import sys
import urllib.request
from pathlib import Path

import pandas as pd

PIPELINE_ROOT = Path(__file__).resolve().parents[1]
if str(PIPELINE_ROOT) not in sys.path:
    sys.path.insert(0, str(PIPELINE_ROOT))

from config import (
    GENCODE_GTF_NAME, GENCODE_TRANSCRIPTS_NAME, MM10_FASTA_NAME, REFERENCE_DIR,
    SACCER3_FASTA_NAME,
)


URLS = {
    f"{GENCODE_GTF_NAME}.gz": "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/gencode.vM25.annotation.gtf.gz",
    f"{GENCODE_TRANSCRIPTS_NAME}.gz": "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/gencode.vM25.transcripts.fa.gz",
    f"{MM10_FASTA_NAME}.gz": "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/GRCm38.primary_assembly.genome.fa.gz",
    f"{SACCER3_FASTA_NAME}.gz": "https://hgdownload.soe.ucsc.edu/goldenPath/sacCer3/bigZips/sacCer3.fa.gz",
}


TOOL_CANDIDATES = {
    "bowtie2-build": ("bowtie2-build", "/opt/anaconda3/envs/cutrun_env/bin/bowtie2-build"),
    "salmon": ("salmon", "/opt/homebrew/bin/salmon", "/opt/anaconda3/bin/salmon"),
    "samtools": ("samtools", "/opt/anaconda3/envs/cutrun_env/bin/samtools"),
}


def require_tool(name: str) -> str:
    """Resolve a documented tool without relying on a caller's PATH."""
    for candidate in TOOL_CANDIDATES[name]:
        value = shutil.which(candidate)
        if value:
            return value
        if Path(candidate).is_file():
            return candidate
    raise FileNotFoundError(f"Required executable not found: {name}")


def download(url: str, destination: Path, dry_run: bool) -> None:
    if destination.is_file():
        return
    print("[DOWNLOAD]", url, "->", destination)
    if not dry_run:
        urllib.request.urlretrieve(url, destination)


def decompress(source: Path, destination: Path, dry_run: bool) -> None:
    if destination.is_file():
        return
    print("[DECOMPRESS]", source, "->", destination)
    if not dry_run:
        with gzip.open(source, "rb") as reader, destination.open("wb") as writer:
            shutil.copyfileobj(reader, writer)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--threads", type=int, default=16)
    parser.add_argument("--full-index", action="store_true", help="Also download FASTAs and build Salmon/Bowtie2 indexes for raw mode")
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args(); REFERENCE_DIR.mkdir(parents=True, exist_ok=True)
    required_urls = URLS if args.full_index else {f"{GENCODE_GTF_NAME}.gz": URLS[f"{GENCODE_GTF_NAME}.gz"]}
    for gz_name, url in required_urls.items():
        gz_path = REFERENCE_DIR / gz_name
        plain_path = REFERENCE_DIR / gz_name.removesuffix(".gz")
        if plain_path.is_file():
            continue
        download(url, gz_path, args.dry_run)
        decompress(gz_path, plain_path, args.dry_run)
    if args.dry_run:
        return
    gtf = REFERENCE_DIR / GENCODE_GTF_NAME; tx2gene = REFERENCE_DIR / "tx2gene.tsv"
    if not tx2gene.is_file():
        rows = []
        with gtf.open() as handle:
            for line in handle:
                if line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) != 9 or fields[2] != "transcript":
                    continue
                attrs = {item.split(" ", 1)[0]: item.split('"')[1] for item in fields[8].split("; ") if '"' in item}
                if "transcript_id" in attrs and "gene_id" in attrs:
                    rows.append((attrs["transcript_id"].split(".")[0], attrs["gene_id"].split(".")[0]))
        pd.DataFrame(sorted(set(rows)), columns=["transcript_id", "gene_id"]).to_csv(tx2gene, sep="\t", index=False)
    if not args.full_index:
        print("[DONE] GENCODE M25 annotation and tx2gene table are available for figure rendering.")
        return
    bowtie2_build = require_tool("bowtie2-build"); salmon = require_tool("salmon"); samtools = require_tool("samtools")
    mouse = REFERENCE_DIR / MM10_FASTA_NAME; yeast = REFERENCE_DIR / SACCER3_FASTA_NAME
    if not (REFERENCE_DIR / "bowtie2_mm10" / "mm10.1.bt2").is_file():
        (REFERENCE_DIR / "bowtie2_mm10").mkdir(exist_ok=True)
        subprocess.run([bowtie2_build, str(mouse), str(REFERENCE_DIR / "bowtie2_mm10" / "mm10")], check=True)
    if not (REFERENCE_DIR / "bowtie2_sacCer3" / "sacCer3.1.bt2").is_file():
        (REFERENCE_DIR / "bowtie2_sacCer3").mkdir(exist_ok=True)
        subprocess.run([bowtie2_build, str(yeast), str(REFERENCE_DIR / "bowtie2_sacCer3" / "sacCer3")], check=True)
    if not (REFERENCE_DIR / "mm10.chrom.sizes").is_file():
        subprocess.run([samtools, "faidx", str(mouse)], check=True)
        with (REFERENCE_DIR / "mm10.chrom.sizes").open("w") as handle:
            for line in (REFERENCE_DIR / f"{MM10_FASTA_NAME}.fai").read_text().splitlines():
                handle.write("\t".join(line.split("\t")[:2]) + "\n")
    if not (REFERENCE_DIR / "salmon_index_M25").is_dir():
        subprocess.run([salmon, "index", "-t", str(REFERENCE_DIR / GENCODE_TRANSCRIPTS_NAME), "-i", str(REFERENCE_DIR / "salmon_index_M25"), "-k", "31", "-p", str(args.threads)], check=True)


if __name__ == "__main__":
    main()
