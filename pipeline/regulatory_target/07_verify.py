#!/usr/bin/env python3
"""Verify that the final workflow emitted every registered publication figure."""

from __future__ import annotations

import csv
import argparse
import re
import sys
from pathlib import Path

PIPELINE_ROOT = Path(__file__).resolve().parents[1]
if str(PIPELINE_ROOT) not in sys.path:
    sys.path.insert(0, str(PIPELINE_ROOT))

from config import CUTRUN_ROOT, REGULATORY_ROOT, RNA_ROOT, RUN_ROOT


EXPECTED_DEG_COUNTS = {
    "MCM3": {"tested": 13196, "significant": 5756, "up": 2656, "down": 3100},
    "NONO": {"tested": 13037, "significant": 653, "up": 446, "down": 207},
    "PSPC1": {"tested": 13037, "significant": 546, "up": 334, "down": 212},
}

def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_gene_set(path: Path, column: str = "gene") -> set[str]:
    return {row[column].strip() for row in read_tsv(path) if row.get(column, "").strip()}


def gencode_symbols() -> dict[str, str]:
    result: dict[str, str] = {}
    gtf = RUN_ROOT / "reference" / "gencode.vM25.annotation.gtf"
    with gtf.open() as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] != "gene":
                continue
            gene_id = re.search(r'gene_id "([^".]+)', fields[8])
            symbol = re.search(r'gene_name "([^"]+)', fields[8])
            if gene_id and symbol:
                result[gene_id.group(1)] = symbol.group(1)
    return result


def verify_rnaseq_and_regulatory(mismatches: list[str]) -> None:
    """Reject stale DEG or regulatory outputs even when every file exists."""
    symbol_map = gencode_symbols()
    direct_sets: dict[str, set[str]] = {}
    for factor, expected in EXPECTED_DEG_COUNTS.items():
        all_rows = read_tsv(RNA_ROOT / "data" / f"Volcano_{factor}_data.tsv")
        significant = [
            row for row in all_rows
            if float(row["adj.P.Val"]) <= 0.05 and abs(float(row["logFC"])) >= 0.28
        ]
        up = {row["gene_id"] for row in significant if float(row["logFC"]) > 0}
        down = {row["gene_id"] for row in significant if float(row["logFC"]) < 0}
        observed = {
            "tested": len(all_rows), "significant": len(significant),
            "up": len(up), "down": len(down),
        }
        if observed != expected:
            mismatches.append(f"{factor} DEG counts: expected {expected}, observed {observed}")
        for path, calculated in (
            (RNA_ROOT / "data" / f"DEGs_{factor}.tsv", {row["gene_id"] for row in significant}),
            (RNA_ROOT / "data" / f"VennDiagram_UP_{factor}_genes.tsv", up),
            (RNA_ROOT / "data" / f"VennDiagram_DOWN_{factor}_genes.tsv", down),
        ):
            written = read_gene_set(path, "gene_id")
            if written != calculated:
                mismatches.append(f"{path.name} is stale or inconsistent")

        deg_symbols = {
            row.get("gene_symbol", "").strip()
            if row.get("gene_symbol", "").strip() not in {"", "NA"}
            else symbol_map.get(row["gene_id"], row["gene_id"])
            for row in significant
        }
        promoter = read_gene_set(CUTRUN_ROOT / "05_promoters" / f"Venn_PromoterGenes_{factor}.tsv")
        expected_direct = promoter & deg_symbols
        direct_sets[factor] = read_gene_set(REGULATORY_ROOT / "data" / f"Venn_target_{factor}_genes.tsv")
        if direct_sets[factor] != expected_direct:
            mismatches.append(f"Venn_target_{factor}_genes.tsv does not use the current {factor} DEG set")

    expected_groups = {
        "A": set.intersection(*direct_sets.values()),
        "B": (direct_sets["MCM3"] & direct_sets["PSPC1"]) - direct_sets["NONO"],
        "C": (direct_sets["NONO"] & direct_sets["PSPC1"]) - direct_sets["MCM3"],
        "D": (direct_sets["MCM3"] & direct_sets["NONO"]) - direct_sets["PSPC1"],
    }
    for group, expected in expected_groups.items():
        observed = read_gene_set(REGULATORY_ROOT / "data" / f"Group{group}_genes.tsv")
        if observed != expected:
            mismatches.append(f"Group{group}.tsv is inconsistent with current regulatory targets")


def clean_resource_forks() -> None:
    """Remove macOS metadata sidecars from generated, reviewer-facing outputs."""
    for root in (RNA_ROOT, CUTRUN_ROOT, REGULATORY_ROOT):
        for path in root.rglob("._*"):
            if path.is_file():
                path.unlink()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()
    if args.dry_run:
        print("[DRY-RUN] Would verify the final figure and data inventories.")
        return
    clean_resource_forks()
    manifest = RUN_ROOT / "docs" / "output_manifest.tsv"
    missing: list[str] = []
    with manifest.open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            path = RUN_ROOT / row["final_relative_output"]
            if not path.is_file() or path.stat().st_size == 0:
                missing.append(str(path.relative_to(RUN_ROOT)))
    required_data = [
        RNA_ROOT / "metadata" / "Samples.tsv", RNA_ROOT / "metadata" / "Salmon_quantifications.tsv",
        RNA_ROOT / "salmon" / "Quantifications_all.tsv",
        RNA_ROOT / "data" / "Volcano_MCM3_data.tsv", RNA_ROOT / "data" / "Volcano_NONO_data.tsv", RNA_ROOT / "data" / "Volcano_PSPC1_data.tsv",
        RNA_ROOT / "data" / "DEGs_MCM3_summary.tsv", RNA_ROOT / "data" / "DEGs_NONO_summary.tsv", RNA_ROOT / "data" / "DEGs_PSPC1_summary.tsv",
        CUTRUN_ROOT / "03_bigwig" / "YeastNormalization.tsv", CUTRUN_ROOT / "data" / "AnalysisParameters.tsv",
        CUTRUN_ROOT / "data" / "PeakLoci" / "Venn_Peaks_loci.tsv", CUTRUN_ROOT / "data" / "Venn_Peaks_counts.tsv",
        CUTRUN_ROOT / "data" / "Pie_PeakDistribution_counts.tsv", CUTRUN_ROOT / "data" / "Venn_Peaks_universe.tsv", CUTRUN_ROOT / "data" / "Peak_profiles_sets.tsv",
        CUTRUN_ROOT / "data" / "Peak_profiles_promoter_overlaps.tsv",
        CUTRUN_ROOT / "data" / "Peak_profiles_promoter_assignments.tsv",
        CUTRUN_ROOT / "05_promoters" / "Venn_PromoterGenes_MCM3.tsv", CUTRUN_ROOT / "05_promoters" / "Venn_PromoterGenes_NONO.tsv", CUTRUN_ROOT / "05_promoters" / "Venn_PromoterGenes_PSPC1.tsv",
        REGULATORY_ROOT / "data" / "promoter_vs_DEG_direction_summary.tsv",
        REGULATORY_ROOT / "data" / "Venn_target_MCM3_genes.tsv", REGULATORY_ROOT / "data" / "Venn_target_NONO_genes.tsv", REGULATORY_ROOT / "data" / "Venn_target_PSPC1_genes.tsv",
        REGULATORY_ROOT / "data" / "Venn_target_shared_genes.tsv",
    ]
    missing.extend(str(path.relative_to(RUN_ROOT)) for path in required_data if not path.is_file() or path.stat().st_size == 0)
    required_dirs = [RNA_ROOT / "salmon", CUTRUN_ROOT / "04_peaks" / "2019", CUTRUN_ROOT / "04_peaks" / "2020"]
    missing.extend(str(path.relative_to(RUN_ROOT)) for path in required_dirs if not path.is_dir())
    if missing:
        raise SystemExit("Missing final outputs:\n" + "\n".join(missing))
    mismatches: list[str] = []
    verify_rnaseq_and_regulatory(mismatches)
    if mismatches:
        raise SystemExit("Inconsistent final outputs:\n" + "\n".join(mismatches))
    print("[PASS] All registered figures are present and RNA-seq/regulatory gene sets are internally consistent.")


if __name__ == "__main__":
    main()
