#!/usr/bin/env python3
"""Validate packaged intermediates and write the local Salmon manifest."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

PIPELINE_ROOT = Path(__file__).resolve().parents[1]
if str(PIPELINE_ROOT) not in sys.path:
    sys.path.insert(0, str(PIPELINE_ROOT))

from config import CUTRUN_ROOT, FACTORS, GENCODE_GTF_NAME, REFERENCE_DIR, RNA_ROOT, RUN_ROOT


def required_intermediates() -> list[Path]:
    paths = [
        REFERENCE_DIR / GENCODE_GTF_NAME,
        RNA_ROOT / "salmon" / "Quantifications_all.tsv",
        CUTRUN_ROOT / "03_bigwig" / "YeastNormalization.tsv",
        CUTRUN_ROOT / "05_consensus" / "NT" / "IgG" / "union.merge.bed",
        CUTRUN_ROOT / "data" / "figure_inputs" / "Venn_Peaks_segments.tsv",
    ]
    for factor in (*FACTORS, "IgG"):
        for replicate in (1, 2):
            library = f"{factor}-{replicate}"
            paths.extend((
                CUTRUN_ROOT / "02_bam" / library / f"{library}.filt.sorted.bam",
                CUTRUN_ROOT / "03_bigwig" / library / f"{library}.bw",
            ))
        paths.extend((
            CUTRUN_ROOT / "03_bigwig" / "mean_intermediates" / f"{factor}_mean.bw",
            CUTRUN_ROOT / "03_bigwig" / "IGV_representation" / f"{factor}_mean.bw",
        ))
    for factor in FACTORS:
        paths.extend((
            CUTRUN_ROOT / "04_peaks" / "pooled" / factor / f"{factor}-pooled_full_peaks.narrowPeak",
            CUTRUN_ROOT / "data" / "figure_inputs" / "promoter_gene_venn" / f"Venn_PromoterGenes_{factor}.tsv",
            CUTRUN_ROOT / "data" / "figure_inputs" / "promoter_gene_venn" / f"Venn_PromoterGenes_{factor}_support.tsv",
        ))
    paths.extend(
        CUTRUN_ROOT / "04_peaks" / f"rep{replicate}" / "IgG" / f"IgG-{replicate}_full_peaks.narrowPeak"
        for replicate in (1, 2)
    )
    return paths


def salmon_relative(path: str) -> Path:
    parts = Path(path).parts
    for index in range(len(parts) - 2):
        if parts[index:index + 2] == ("deg_work", "salmon"):
            return Path(*parts[index + 2:])
    raise ValueError(f"Could not locate deg_work/salmon in {path}")


def write_quant_manifest(dry_run: bool) -> int:
    samples = pd.read_csv(RNA_ROOT / "metadata" / "Samples.tsv", sep="\t", dtype=str)
    selected = samples["selected"].str.lower().eq("true")
    context = samples["legacy_nono_pspc1_context"].str.lower().eq("true")
    samples = samples.loc[selected | context].copy()

    accepted = pd.read_csv(RNA_ROOT / "salmon" / "Quantifications_all.tsv", sep="\t", dtype=str)
    accepted["r1_basename"] = accepted["r1_path"].map(lambda value: Path(value).name)
    samples["r1_basename"] = samples["r1"].map(lambda value: Path(value).name)
    if accepted["r1_basename"].duplicated().any():
        raise RuntimeError("Accepted Salmon metadata contains duplicate R1 basenames")
    merged = samples.merge(
        accepted[["r1_basename", "replicate_id", "quant_sf"]],
        on="r1_basename", how="left", validate="one_to_one",
    )
    if merged["quant_sf"].isna().any():
        missing = merged.loc[merged["quant_sf"].isna(), "r1_basename"].tolist()
        raise RuntimeError(f"Selected FASTQs lack Salmon quantifications: {missing}")
    merged["quant_sf"] = merged["quant_sf"].map(
        lambda value: str(RNA_ROOT / "salmon" / salmon_relative(value))
    )
    missing = [path for path in merged["quant_sf"] if not Path(path).is_file()]
    if missing:
        raise FileNotFoundError(f"Missing packaged quant.sf files: {missing[:3]}")

    output = RNA_ROOT / "metadata" / "Salmon_quantifications.tsv"
    columns = [
        "sample_id", "experiment", "condition", "target", "shRNA", "bioreplicate",
        "r1", "r2", "selected", "analysis", "legacy_nono_pspc1_context",
        "legacy_context_order", "replicate_id", "quant_sf",
    ]
    if not dry_run:
        merged[columns].to_csv(output, sep="\t", index=False)
    return len(merged)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()

    missing = [path for path in required_intermediates() if not path.is_file()]
    if missing:
        raise FileNotFoundError(
            "Packaged accepted intermediates are incomplete. Use --source raw to rebuild them. "
            f"First missing path: {missing[0]}"
        )
    count = write_quant_manifest(args.dry_run)
    if not args.dry_run:
        output = RUN_ROOT / "docs" / "accepted_intermediate_manifest.tsv"
        pd.DataFrame([
            {"intermediate": "RNA-seq Salmon quantifications", "location": "deg_work/salmon"},
            {"intermediate": "2020 CUT&RUN filtered alignments", "location": "cutrun_work/02_bam"},
            {"intermediate": "2020 yeast-normalized tracks", "location": "cutrun_work/03_bigwig"},
            {"intermediate": "2020 matched-IgG peak calls", "location": "cutrun_work/04_peaks"},
            {"intermediate": "2020 peak and promoter membership", "location": "cutrun_work/data/figure_inputs"},
        ]).to_csv(output, sep="\t", index=False)
    status = "DRY-RUN" if args.dry_run else "DONE"
    print(f"[{status}] Validated packaged intermediates; Salmon manifest contains {count} samples.")


if __name__ == "__main__":
    main()
