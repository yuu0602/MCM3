#!/usr/bin/env python3
"""Render RNA-seq Venn diagrams."""

from __future__ import annotations

import argparse
import os
import subprocess
import sys
from pathlib import Path

import pandas as pd

PIPELINE_ROOT = Path(__file__).resolve().parents[1]
if str(PIPELINE_ROOT) not in sys.path:
    sys.path.insert(0, str(PIPELINE_ROOT))

from config import FACTORS, PIPELINE_ROOT, RNA_FDR_MAX, RNA_LOGFC_MIN, RNA_ROOT


HELPER = Path(__file__).resolve().parent / "figures_rendering" / "venn.py"
OUTPUT = RNA_ROOT / "visuals"


def genes_for_direction(factor: str, direction: str) -> set[str]:
    table = pd.read_csv(RNA_ROOT / "data" / f"Volcano_{factor}_data.tsv", sep="\t", dtype={"gene_id": str})
    keep = (pd.to_numeric(table["adj.P.Val"], errors="coerce") <= RNA_FDR_MAX) & (
        pd.to_numeric(table["logFC"], errors="coerce").abs() >= RNA_LOGFC_MIN
    )
    logfc = pd.to_numeric(table["logFC"], errors="coerce")
    keep &= logfc.gt(0) if direction == "UP" else logfc.lt(0)
    return set(table.loc[keep, "gene_id"].dropna().astype(str))


def regions(sets: dict[str, set[str]]) -> dict[str, int]:
    a, b, c = (sets[factor] for factor in FACTORS)
    return {
        "n100": len(a - b - c), "n010": len(b - a - c), "n001": len(c - a - b),
        "n110": len((a & b) - c), "n101": len((a & c) - b),
        "n011": len((b & c) - a), "n111": len(a & b & c),
    }


def render_direction(direction: str) -> None:
    if not HELPER.is_file():
        raise FileNotFoundError(HELPER)
    sets = {factor: genes_for_direction(factor, direction) for factor in FACTORS}
    counts = regions(sets)
    title = "Unregulated genes" if direction == "UP" else "Downregulated genes"
    OUTPUT.mkdir(parents=True, exist_ok=True)
    environment = os.environ.copy()
    environment.setdefault("MPLCONFIGDIR", str(RNA_ROOT.parent / ".cache" / "matplotlib"))
    for suffix, hide_numbers in (("", False), ("_noNumbers", True)):
        output = OUTPUT / f"VennDiagram_{direction}{suffix}.png"
        command = [
            sys.executable, str(HELPER), "--out", str(output), "--title", title,
            "--a-name", "MCM3", "--b-name", "NONO", "--c-name", "PSPC1",
            "--a-total", str(len(sets["MCM3"])), "--b-total", str(len(sets["NONO"])), "--c-total", str(len(sets["PSPC1"])),
            *[item for pair in counts.items() for item in (f"--{pair[0]}", str(pair[1]))],
        ]
        if hide_numbers:
            command.append("--hide-numbers")
        subprocess.run(command, check=True, env=environment)
    pd.DataFrame([
        {"factor": factor, "direction": direction, "gene_id": gene}
        for factor, genes in sets.items() for gene in sorted(genes)
    ]).to_csv(RNA_ROOT / "data" / f"VennDiagram_{direction}_data.tsv", sep="\t", index=False)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()
    if args.dry_run:
        print("[DRY-RUN] Would render RNA-seq UP and DOWN Venn diagrams.")
        return
    for direction in ("UP", "DOWN"):
        render_direction(direction)
    print("[DONE] RNA-seq Venn diagrams rendered.")


if __name__ == "__main__":
    main()
