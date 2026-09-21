"""GENCODE and deepTools matrix helpers used by CUT&RUN figures."""

import gzip
import re
from pathlib import Path

import numpy as np
import pandas as pd


def gene_table(gtf: Path) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    with gtf.open() as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] != "gene":
                continue
            gene_id = re.search(r'gene_id "([^"]+)"', fields[8])
            gene_name = re.search(r'gene_name "([^"]+)"', fields[8])
            if gene_id:
                rows.append({
                    "chrom": fields[0],
                    "start": int(fields[3]) - 1,
                    "end": int(fields[4]),
                    "strand": fields[6],
                    "gene_id": gene_id.group(1),
                    "gene_name": gene_name.group(1) if gene_name else gene_id.group(1),
                })
    table = pd.DataFrame(rows)
    if table.empty:
        raise RuntimeError(f"No gene records found in {gtf}")
    table["gene_length"] = table["end"] - table["start"]
    table["canonical_chromosome"] = table["chrom"].str.match(r"^chr(?:[1-9]|1[0-9]|X|Y|M)$")
    return (
        table.sort_values(
            ["gene_name", "canonical_chromosome", "gene_length", "gene_id"],
            ascending=[True, False, False, True],
        )
        .drop_duplicates("gene_name", keep="first")
        .copy()
    )


def write_gene_bed(gtf: Path, symbols: set[str], output: Path) -> int:
    table = gene_table(gtf)
    selected = table.loc[table.gene_name.isin(symbols)].copy()
    missing = symbols - set(selected.gene_name)
    if missing:
        raise RuntimeError(f"{len(missing)} promoter genes are absent from GENCODE M25")
    selected.insert(5, "score", 0)
    selected = selected[["chrom", "start", "end", "gene_id", "score", "strand"]]
    selected = selected.sort_values(["chrom", "start", "end", "gene_id"])
    output.parent.mkdir(parents=True, exist_ok=True)
    selected.to_csv(output, sep="\t", header=False, index=False)
    return len(selected)


def read_matrix(path: Path, n_tracks: int) -> tuple[list[str], np.ndarray]:
    names: list[str] = []
    rows: list[np.ndarray] = []
    with gzip.open(path, "rt") as handle:
        for line in handle:
            if not line.strip() or line.startswith(("@", "#")):
                continue
            fields = line.rstrip("\n").split("\t")
            values = np.array([0.0 if value == "nan" else float(value) for value in fields[6:]])
            if len(values) % n_tracks:
                raise RuntimeError(f"Unexpected matrix width in {path}")
            names.append(fields[3])
            rows.append(values.reshape(n_tracks, len(values) // n_tracks))
    if not rows:
        raise RuntimeError(f"No regions found in {path}")
    return names, np.stack(rows)
