#!/usr/bin/env python3
"""Verify that the final workflow emitted every registered publication figure."""

from __future__ import annotations

import csv
import argparse
import gzip
import re
import sys
from pathlib import Path

PIPELINE_ROOT = Path(__file__).resolve().parents[1]
if str(PIPELINE_ROOT) not in sys.path:
    sys.path.insert(0, str(PIPELINE_ROOT))

from config import CUTRUN_ROOT, REGULATORY_ROOT, RNA_ROOT, RUN_ROOT


EXPECTED_DEG_COUNTS = {
    "MCM3": {"tested": 13197, "significant": 5757, "up": 2657, "down": 3100},
    "NONO": {"tested": 13038, "significant": 653, "up": 446, "down": 207},
    "PSPC1": {"tested": 13038, "significant": 549, "up": 333, "down": 216},
}

EXPECTED_CUTRUN_COUNTS = {
    "MCM3": {"peak_segments": 11184, "promoter_bound_genes": 2945},
    "NONO": {"peak_segments": 2500, "promoter_bound_genes": 745},
    "PSPC1": {"peak_segments": 20666, "promoter_bound_genes": 4903},
}

EXPECTED_PROMOTER_REGIONS = {
    "100_only_MCM3": 248,
    "010_only_NONO": 36,
    "001_only_PSPC1": 2128,
    "110_MCM3_NONO": 12,
    "101_MCM3_PSPC1": 2078,
    "011_NONO_PSPC1": 90,
    "111_MCM3_NONO_PSPC1": 607,
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


def verify_cutrun(mismatches: list[str]) -> None:
    support = {
        row["factor"]: {
            "peak_segments": int(row["peak_segments"]),
            "promoter_bound_genes": int(row["promoter_bound_genes"]),
        }
        for row in read_tsv(CUTRUN_ROOT / "data" / "Venn_Peaks_support.tsv")
    }
    if support != EXPECTED_CUTRUN_COUNTS:
        mismatches.append(
            f"CUT&RUN counts: expected {EXPECTED_CUTRUN_COUNTS}, observed {support}"
        )
    regions: dict[str, int] = {}
    path = (
        CUTRUN_ROOT
        / "data"
        / "figure_inputs"
        / "promoter_gene_venn"
        / "Venn_PromoterGenes_regions.tsv"
    )
    for row in read_tsv(path):
        regions[row["region"]] = regions.get(row["region"], 0) + 1
    if regions != EXPECTED_PROMOTER_REGIONS:
        mismatches.append(
            f"promoter-region counts: expected {EXPECTED_PROMOTER_REGIONS}, observed {regions}"
        )
    verify_promoter_profile_groups(mismatches)


def verify_promoter_profile_groups(mismatches: list[str]) -> None:
    factors = ("MCM3", "NONO", "PSPC1")
    sets = {factor: read_gene_set(CUTRUN_ROOT / "05_promoters" / f"Venn_PromoterGenes_{factor}.tsv") for factor in factors}
    data = CUTRUN_ROOT / "data"
    rows = read_tsv(data / "MCM3_Peak_profiles_genes.tsv")
    listed = [row["gene"] for row in rows]
    if len(listed) != len(set(listed)) or set(listed) != sets["MCM3"]:
        mismatches.append("MCM3 gene-profile groups must partition the final MCM3 promoter-bound set")
    for row in rows:
        expected = sum(1 << i for i, factor in enumerate(factors) if row["gene"] in sets[factor])
        if int(row["promoter_binding_mask"]) != expected:
            mismatches.append(f"Incorrect gene-level binding group: {row['gene']}")
    triple_rows = read_tsv(data / "MCM3_NONO_PSPC1_peak_profiles_genes.tsv")
    triple = [row["gene"] for row in triple_rows]
    expected_triple = set.intersection(*sets.values())
    if len(triple) != len(set(triple)) or set(triple) != expected_triple:
        mismatches.append("Triple profile gene list disagrees with the promoter-bound intersection")
    matrices = {
        f"MCM3_Peak_profiles_region_{mask:03b}_matrix.gz": {
            row["gene"] for row in rows if int(row["promoter_binding_mask"]) == mask
        } for mask in (1, 3, 5, 7)
    }
    matrices["MCM3_NONO_PSPC1_peak_profiles_matrix.gz"] = expected_triple
    for name, expected in matrices.items():
        path = data / "figure_inputs" / "additional_panels" / name
        if not path.is_file():
            mismatches.append(f"Missing profile matrix: {name}")
            continue
        with gzip.open(path, "rt") as handle:
            genes = [line.rstrip("\n").split("\t")[3] for line in handle if line.strip() and not line.startswith(("@", "#"))]
        if len(genes) != len(set(genes)) or set(genes) != expected:
            mismatches.append(f"{name}: expected one matrix row per promoter-bound gene")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()
    if args.dry_run:
        print("[DRY-RUN] Would verify the final figure and data inventories.")
        return
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
        CUTRUN_ROOT / "data" / "MCM3_Peak_profiles_genes.tsv",
        CUTRUN_ROOT / "data" / "MCM3_NONO_PSPC1_peak_profiles_genes.tsv",
        CUTRUN_ROOT / "05_promoters" / "Venn_PromoterGenes_MCM3.tsv", CUTRUN_ROOT / "05_promoters" / "Venn_PromoterGenes_NONO.tsv", CUTRUN_ROOT / "05_promoters" / "Venn_PromoterGenes_PSPC1.tsv",
        REGULATORY_ROOT / "data" / "promoter_vs_DEG_direction_summary.tsv",
        REGULATORY_ROOT / "data" / "Venn_target_MCM3_genes.tsv", REGULATORY_ROOT / "data" / "Venn_target_NONO_genes.tsv", REGULATORY_ROOT / "data" / "Venn_target_PSPC1_genes.tsv",
        REGULATORY_ROOT / "data" / "Venn_target_shared_genes.tsv",
    ]
    missing.extend(str(path.relative_to(RUN_ROOT)) for path in required_data if not path.is_file() or path.stat().st_size == 0)
    required_dirs = [RNA_ROOT / "salmon", CUTRUN_ROOT / "04_peaks"]
    missing.extend(str(path.relative_to(RUN_ROOT)) for path in required_dirs if not path.is_dir())
    if missing:
        raise SystemExit("Missing final outputs:\n" + "\n".join(missing))
    mismatches: list[str] = []
    verify_cutrun(mismatches)
    verify_rnaseq_and_regulatory(mismatches)
    if mismatches:
        raise SystemExit("Inconsistent final outputs:\n" + "\n".join(mismatches))
    print("[PASS] All registered figures are present and RNA-seq/regulatory gene sets are internally consistent.")


if __name__ == "__main__":
    main()
