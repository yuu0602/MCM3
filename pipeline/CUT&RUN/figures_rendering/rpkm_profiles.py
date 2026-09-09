#!/usr/bin/env python3
"""Render promoter-binding RPKM boxplots for each CUT&RUN factor."""

from __future__ import annotations

import csv
import gzip
import math
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Iterable


CUTRUN_ROOT = Path.cwd() / "cutrun_work_20191012"
FIGURE_DATA = CUTRUN_ROOT / "data" / "figure_inputs"
OUTDIR = FIGURE_DATA / "rpkm_profiles" / "intermediate"
TMPDIR = OUTDIR / "regions"
MATRIX_TMPDIR = Path("/tmp/mcm3_cutrun_rpkm_matrices")
REGION_TABLE = FIGURE_DATA / "promoter_gene_venn" / "Venn_PromoterGenes_regions.tsv"
GENE_BED = CUTRUN_ROOT / "data" / "GeneBodies_M25.bed6"
COMPUTE_MATRIX = Path("/opt/anaconda3/envs/cutrun_env/bin/computeMatrix")
LOCAL_BIGWIG_CACHE = Path("/tmp/mcm3_cutrun_bigwig_cache_20191012")
LOCAL_GENE_BED = LOCAL_BIGWIG_CACHE / "genes.norm.sorted.bed6"

OUTPUT = OUTDIR / "MCM3_RPKM.png"
OUTPUT_NONO = OUTDIR / "NONO_RPKM.png"
OUTPUT_PSPC1 = OUTDIR / "PSPC1_RPKM.png"
SUMMARY = OUTDIR / "RPKM_summary.tsv"

ASSAYS = ("MCM3", "NONO", "PSPC1")
BIGWIGS = {
    "MCM3": CUTRUN_ROOT / "03_bigwig/NT_rep1_MCM3/NT_rep1_MCM3.bw",
    "NONO": CUTRUN_ROOT / "03_bigwig/NT_rep1_NONO/NT_rep1_NONO.bw",
    "PSPC1": CUTRUN_ROOT / "03_bigwig/NT_rep1_PSPC1/NT_rep1_PSPC1.bw",
}
REGION_SPECS = (
    ("111", "111_MCM3_NONO_PSPC1", "MCM3-NONO-PSPC1 co-binding", 3473),
    ("110", "110_MCM3_NONO", "MCM3-NONO co-binding", 177),
    ("101", "101_MCM3_PSPC1", "MCM3-PSPC1 co-binding", 2307),
    ("100", "100_only_MCM3", "MCM3-specific binding", 1113),
)
REGION_STEMS = {
    "111": "RPKM_MCM3_NONO_PSPC1",
    "110": "RPKM_MCM3_NONO",
    "101": "RPKM_MCM3_PSPC1",
    "100": "RPKM_MCM3_only",
}
FACTOR_SPECS = (
    ("MCM3", 0, OUTPUT),
    ("NONO", 1, OUTPUT_NONO),
    ("PSPC1", 2, OUTPUT_PSPC1),
)

UPSTREAM_BP = 3000
BODY_BP = 2000
DOWNSTREAM_BP = 1000
BIN_SIZE = 25
TSS_BINS = 4
EXPECTED_BINS = (UPSTREAM_BP + BODY_BP + DOWNSTREAM_BP) // BIN_SIZE
DISPLAY_SCALE = 1000.0 / BIN_SIZE
SIGNAL_LABEL = "Log2(RPKM)"
Y_MAX = 6


def die(msg: str) -> None:
    raise SystemExit(f"[ERROR] {msg}")


def mean(values: Iterable[float]) -> float:
    vals = list(values)
    return sum(vals) / float(len(vals)) if vals else 0.0


def count_matrix_rows(matrix_path: Path) -> int:
    with gzip.open(matrix_path, "rt") as handle:
        return sum(1 for line in handle if line.strip() and not line.startswith(("@", "#")))


def count_bed_rows(bed_path: Path) -> int:
    with bed_path.open() as handle:
        return sum(1 for line in handle if line.strip())


def matrix_is_valid(matrix_path: Path, expected_rows: int) -> bool:
    if not matrix_path.is_file() or matrix_path.stat().st_size == 0:
        return False
    try:
        with gzip.open(matrix_path, "rt") as handle:
            for line in handle:
                if not line.strip() or line.startswith(("@", "#")):
                    continue
                signal_columns = len(line.rstrip("\n").split("\t")) - 6
                return signal_columns == len(ASSAYS) * EXPECTED_BINS and count_matrix_rows(matrix_path) == expected_rows
    except OSError:
        return False
    return False


def resolve_bigwigs() -> dict[str, Path]:
    cached = {assay: LOCAL_BIGWIG_CACHE / source.name for assay, source in BIGWIGS.items()}
    if all(path.is_file() and path.stat().st_size == BIGWIGS[assay].stat().st_size for assay, path in cached.items()):
        print(f"[OK] Using byte-matched local bigWig cache: {LOCAL_BIGWIG_CACHE}")
        return cached
    missing = [str(path) for path in BIGWIGS.values() if not path.is_file() or path.stat().st_size == 0]
    if missing:
        die("Missing CUT&RUN bigWig(s):\n" + "\n".join(missing))
    print("[OK] Using source bigWigs")
    return BIGWIGS


def resolve_gene_bed() -> Path:
    if LOCAL_GENE_BED.is_file() and GENE_BED.is_file() and LOCAL_GENE_BED.stat().st_size == GENE_BED.stat().st_size:
        return LOCAL_GENE_BED
    if not GENE_BED.is_file() or GENE_BED.stat().st_size == 0:
        die(f"Missing normalized gene-body BED: {GENE_BED}")
    return GENE_BED


def load_region_genes() -> dict[str, set[str]]:
    if not REGION_TABLE.is_file() or REGION_TABLE.stat().st_size == 0:
        die(f"Missing promoter-gene Venn-region table: {REGION_TABLE}")
    result: dict[str, set[str]] = {}
    with REGION_TABLE.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if set(reader.fieldnames or []) != {"region", "gene"}:
            die(f"Unexpected region-table columns: {reader.fieldnames}")
        wanted = {table_region for _, table_region, _, _ in REGION_SPECS}
        for row in reader:
            table_region = row["region"]
            gene = row["gene"]
            if table_region in wanted and gene:
                result.setdefault(table_region, set()).add(gene)

    for _, table_region, _, expected_n in REGION_SPECS:
        n_genes = len(result.get(table_region, set()))
        if n_genes != expected_n:
            die(f"{table_region} contains {n_genes:,} genes; expected {expected_n:,}.")
    return result


def write_region_beds(region_genes: dict[str, set[str]], annotation_bed: Path) -> dict[str, Path]:
    records: dict[str, list[str]] = {}
    with annotation_bed.open() as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 6:
                records.setdefault(fields[3], []).append(line)

    TMPDIR.mkdir(parents=True, exist_ok=True)
    outputs: dict[str, Path] = {}
    for code, table_region, _, expected_n in REGION_SPECS:
        genes = region_genes[table_region]
        missing = sorted(gene for gene in genes if gene not in records)
        if missing:
            die(f"{code} contains {len(missing):,} genes absent from the GTF BED, including: {', '.join(missing[:20])}")
        out_bed = TMPDIR / f"{REGION_STEMS[code]}_regions.bed"
        rows = [line for gene in genes for line in records[gene]]
        rows.sort(key=lambda line: (line.split("\t", 3)[0], int(line.split("\t", 3)[1]), int(line.split("\t", 3)[2]), line.split("\t", 4)[3]))
        with out_bed.open("w") as handle:
            handle.writelines(rows)
        n_unique = len({line.rstrip("\n").split("\t")[3] for line in rows})
        if n_unique != expected_n:
            die(f"{code} gene BED has {n_unique:,} unique genes; expected {expected_n:,}.")
        outputs[code] = out_bed
    return outputs


def compute_matrix(region_bed: Path, output_matrix: Path, bigwigs: dict[str, Path]) -> None:
    command = [
        str(COMPUTE_MATRIX),
        "scale-regions",
        "-S",
        *[str(bigwigs[assay]) for assay in ASSAYS],
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
        "-o",
        str(output_matrix),
    ]
    print("[CMD] " + " ".join(command))
    subprocess.run(command, check=True)


def ensure_matrices(region_beds: dict[str, Path], bigwigs: dict[str, Path]) -> dict[str, Path]:
    MATRIX_TMPDIR.mkdir(parents=True, exist_ok=True)
    outputs: dict[str, Path] = {}
    for code, _, _, _ in REGION_SPECS:
        matrix_path = MATRIX_TMPDIR / f"{REGION_STEMS[code]}_matrix.gz"
        expected_rows = count_bed_rows(region_beds[code])
        if matrix_is_valid(matrix_path, expected_rows):
            print(f"[OK] Reusing validated matrix: {matrix_path}")
        else:
            compute_matrix(region_beds[code], matrix_path, bigwigs)
            if not matrix_is_valid(matrix_path, expected_rows):
                die(f"Invalid matrix after computeMatrix: {matrix_path}")
        outputs[code] = matrix_path
    return outputs


def parse_region_values(matrix_path: Path, assay_idx: int) -> list[float]:
    n_up = UPSTREAM_BP // BIN_SIZE
    out: list[float] = []
    with gzip.open(matrix_path, "rt") as handle:
        for line in handle:
            if not line.strip() or line.startswith(("@", "#")):
                continue
            fields = line.rstrip("\n").split("\t")
            signal = fields[6:]
            assays = [signal[index * EXPECTED_BINS : (index + 1) * EXPECTED_BINS] for index in range(len(ASSAYS))]
            values = [0.0 if token == "nan" else float(token) for token in assays[assay_idx][n_up : n_up + TSS_BINS]]
            out.append(math.log2(mean(values) * DISPLAY_SCALE + 1.0))
    if not out:
        die(f"No usable rows in matrix: {matrix_path}")
    return out


def write_values_tsv(tsv_path: Path, assay_idx: int, matrices: dict[str, Path]) -> None:
    with tsv_path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["group", "value"])
        for code, _, label, _ in REGION_SPECS:
            for value in parse_region_values(matrices[code], assay_idx):
                writer.writerow([label, f"{value:.10f}"])


def write_summary(matrices: dict[str, Path], region_beds: dict[str, Path]) -> None:
    with SUMMARY.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["region_code", "venn_region", "group_label", "n_unique_genes", "n_gene_body_records", "matrix"])
        for code, table_region, label, expected_n in REGION_SPECS:
            writer.writerow([code, table_region, label, expected_n, count_bed_rows(region_beds[code]), matrices[code]])


def r_plot(values_tsv: Path, output_png: Path) -> None:
    r_script = r'''
args <- commandArgs(trailingOnly = TRUE)
in_tsv <- args[1]
out_png <- args[2]
signal_label <- args[3]
y_max <- as.numeric(args[4])

df <- read.delim(in_tsv, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
lev <- c("MCM3-NONO-PSPC1 co-binding",
         "MCM3-NONO co-binding",
         "MCM3-PSPC1 co-binding",
         "MCM3-specific binding")
df$group <- factor(df$group, levels = lev)

cols <- c("#cf111f", "#f2a010", "#4f7db1", "#777777")

png(out_png, width = 2400, height = 1800, res = 300, bg = "white")
par(mar = c(10, 7, 2, 2) + 0.1, cex.axis = 1.25, cex.lab = 1.45, font.axis = 2, font.lab = 2)

boxplot(
  value ~ group, data = df,
  col = cols, border = "#222222", lwd = 2.0,
  ylab = signal_label, xlab = "",
  xaxt = "n", yaxt = "n", outline = FALSE,
  ylim = c(0, y_max),
  whisklty = 1,
  staplelty = 1,
  cex.lab = 1.45, font.lab = 2
)

axis(1, at = seq_along(lev), labels = FALSE, tick = FALSE)
lab <- par("usr")[3] - 0.32
text(
  x = seq_along(lev),
  y = rep(lab, length(lev)),
  labels = lev,
  srt = 35,
  xpd = TRUE,
  adj = 1,
  cex = 0.9,
  font = 2
)

axis(2, at = seq(0, y_max, by = 1), labels = seq(0, y_max, by = 1), las = 1, cex.axis = 1.25, font = 2)
dev.off()
'''
    proc = subprocess.run(
        ["Rscript", "-", str(values_tsv), str(output_png), SIGNAL_LABEL, str(Y_MAX)],
        input=r_script,
        text=True,
        capture_output=True,
        cwd=str(Path.cwd()),
    )
    if proc.returncode != 0:
        sys.stderr.write(proc.stdout or "")
        sys.stderr.write(proc.stderr or "")
        die("R plotting step failed.")


def main() -> None:
    required = [CUTRUN_ROOT, REGION_TABLE, GENE_BED, COMPUTE_MATRIX]
    missing = [str(path) for path in required if not path.exists() or path.stat().st_size == 0]
    if missing:
        die("Missing required input(s):\n" + "\n".join(missing))

    region_genes = load_region_genes()
    region_beds = write_region_beds(region_genes, resolve_gene_bed())
    matrices = ensure_matrices(region_beds, resolve_bigwigs())
    write_summary(matrices, region_beds)

    with tempfile.TemporaryDirectory(prefix="figure4f_") as td:
        td_path = Path(td)
        for factor_name, assay_idx, out_png in FACTOR_SPECS:
            values_tsv = td_path / f"RPKM_values_{factor_name}.tsv"
            write_values_tsv(values_tsv, assay_idx, matrices)
            r_plot(values_tsv, out_png)
            print(f"[OK] Wrote: {out_png}")
    print(f"[OK] Wrote: {SUMMARY}")


if __name__ == "__main__":
    main()
