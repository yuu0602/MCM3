#!/usr/bin/env python3
"""Render all bulk RNA-seq publication figures from packaged quantifications."""

from __future__ import annotations

import argparse
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

from config import FACTORS, REFERENCE_DIR, RNA_FDR_MAX, RNA_LOGFC_MIN, RNA_ROOT


VENN_RENDERER = SCRIPT_DIR / "figures_rendering" / "venn.py"
OUTPUT = RNA_ROOT / "visuals"
PUBLICATION_FIGURES = OUTPUT / "publication_figures"


def executable(name: str) -> str:
    path = shutil.which(name)
    if path:
        return path
    raise FileNotFoundError(name)


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
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--publication-figures", action="store_true", help="Also render text-free PNGs in deg_work/visuals/publication_figures")
    args = parser.parse_args()

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
        return

    subprocess.run(r_command, check=True)
    for direction in ("UP", "DOWN"):
        render_venn(direction, OUTPUT)
    if args.publication_figures:
        for direction in ("UP", "DOWN"):
            render_venn(direction, PUBLICATION_FIGURES, no_text=True)
        print(f"[DONE] Bulk RNA-seq publication figures: {PUBLICATION_FIGURES}")
    print(f"[DONE] Bulk RNA-seq figures: {RNA_ROOT / 'visuals'}")


if __name__ == "__main__":
    main()
