#!/usr/bin/env python3
"""Create final analysis manifests directly from original FASTQs.

L001-L004 are retained as independent biological samples.  This stage never
reads a result file from any earlier reanalysis.
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import pandas as pd

PIPELINE_ROOT = Path(__file__).resolve().parents[1]
if str(PIPELINE_ROOT) not in sys.path:
    sys.path.insert(0, str(PIPELINE_ROOT))

from config import (
    CUTRUN_ROOT,
    RAW_CUTRUN_2019,
    RAW_CUTRUN_2020,
    RAW_RNASETS,
    REFERENCE_DIR,
    RNA_ROOT,
)


def require(path: Path, label: str) -> None:
    if not path.exists():
        raise FileNotFoundError(f"Missing {label}: {path}")


def lane_from_name(name: str) -> str:
    hit = re.search(r"_L(\d{3})_", name)
    if hit is None:
        raise ValueError(f"Cannot infer L00X biological replicate from {name}")
    return f"L{hit.group(1)}"


def classify_rna(experiment: str, parent: str) -> tuple[str, str, str]:
    text = parent.lower()
    if experiment == "20161213":
        if "mcm3kd_nt" in text:
            return "NT", "MCM3", "NT"
        if "sh1" in text:
            return "KD", "MCM3", "sh1"
        if "sh2" in text:
            return "KD", "MCM3", "sh2"
    if experiment == "201709":
        if parent.startswith("NT-201709-"):
            return "NT", "MCM3", "NT"
        if parent.startswith("sh3-6-"):
            return "KD", "MCM3", "sh3_6"
        if parent.startswith("sh3-10-"):
            return "KD", "MCM3", "sh3_10"
    if experiment == "20180718":
        if parent.startswith("NT-24_L"):
            return "NT", "NT", "NT"
        hit = re.match(r"^(Nono|Pspc1)_sh(\d+)-24_L", parent, flags=re.I)
        if hit:
            return "KD", hit.group(1).upper(), f"sh{hit.group(2)}"
    raise ValueError(f"Unrecognized RNA-seq sample directory: {experiment}/{parent}")


def include_rna(experiment: str, condition: str, target: str, shrna: str, lane: str) -> tuple[bool, str]:
    """Encode the final RNA-seq sample selection."""
    if target == "MCM3":
        keep = experiment == "20161213" and (
            (condition == "NT" and lane in {"L001", "L002", "L003", "L004"})
            or (shrna == "sh1" and lane in {"L001", "L002", "L003", "L004"})
            or (shrna == "sh2" and lane in {"L002", "L003"})
        )
        return keep, "MCM3_NT_vs_sh1_sh2" if keep else "not_selected"
    if target in {"NONO", "PSPC1"}:
        if experiment != "20180718":
            return False, "not_selected"
        if condition == "NT":
            keep = lane in {"L001", "L002", "L003", "L004"}
        elif target == "NONO":
            keep = (shrna, lane) in {("sh455", "L001"), ("sh455", "L003"), ("sh777", "L001"), ("sh777", "L003")}
        else:
            keep = (shrna, lane) in {("sh197", "L001"), ("sh197", "L002"), ("sh873", "L002"), ("sh873", "L003")}
        return keep, f"{target}_NT_vs_pooledKD" if keep else "not_selected"
    if target == "NT" and experiment == "20180718":
        return lane in {"L001", "L002", "L003", "L004"}, "shared_20180718_NT"
    return False, "not_selected"


def legacy_nono_pspc1_context(
    experiment: str, condition: str, target: str, shrna: str, lane: str
) -> tuple[bool, int | None]:
    """Return the retained reanalysis11 normalization cohort and its order.

    NONO and PSPC1 were accepted from the original per-target analysis.  That
    analysis calculated TMM factors and its expression filter on this fixed
    22-library cohort before fitting the two 20180718 target contrasts.  The
    MCM3 contrast itself is still fitted independently from the current
    sh1+sh2 sample set.
    """
    key = (experiment, condition, target, shrna, lane)
    order: list[tuple[str, str, str, str, str]] = []
    order.extend(("20161213", "KD", "MCM3", "sh1", value) for value in ("L001", "L002", "L003", "L004"))
    order.extend(("201709", "KD", "MCM3", "sh3_6", value) for value in ("L001", "L003"))
    order.extend(("20161213", "NT", "MCM3", "NT", value) for value in ("L002", "L004"))
    order.extend(("201709", "NT", "MCM3", "NT", value) for value in ("L001", "L002"))
    order.extend(("20180718", "KD", "NONO", "sh455", value) for value in ("L001", "L003"))
    order.extend(("20180718", "KD", "NONO", "sh777", value) for value in ("L001", "L003"))
    order.extend(("20180718", "NT", "NT", "NT", value) for value in ("L001", "L002", "L003", "L004"))
    order.extend(("20180718", "KD", "PSPC1", "sh197", value) for value in ("L001", "L002"))
    order.extend(("20180718", "KD", "PSPC1", "sh873", value) for value in ("L002", "L003"))
    try:
        return True, order.index(key) + 1
    except ValueError:
        return False, None


def build_rna_manifest() -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for experiment, root in RAW_RNASETS.items():
        require(root, f"RNA-seq FASTQ root for {experiment}")
        for r1 in sorted(root.rglob("*_R1_001.fastq.gz")):
            if "Undetermined" in str(r1) or "Unindexed_Reads" in str(r1) or any(part.startswith("._") for part in r1.parts):
                continue
            r2 = r1.with_name(r1.name.replace("_R1_001.fastq.gz", "_R2_001.fastq.gz"))
            if not r2.is_file():
                raise FileNotFoundError(f"Missing paired read for {r1}")
            condition, target, shrna = classify_rna(experiment, r1.parent.name)
            lane = lane_from_name(r1.name)
            selected, analysis = include_rna(experiment, condition, target, shrna, lane)
            legacy_context, legacy_order = legacy_nono_pspc1_context(
                experiment, condition, target, shrna, lane
            )
            rows.append(
                {
                    "sample_id": f"{experiment}_{r1.parent.name}_{lane}",
                    "experiment": experiment,
                    "condition": condition,
                    "target": target,
                    "shRNA": shrna,
                    "bioreplicate": lane,
                    "r1": str(r1.resolve()),
                    "r2": str(r2.resolve()),
                    "selected": selected,
                    "analysis": analysis,
                    "legacy_nono_pspc1_context": legacy_context,
                    "legacy_context_order": legacy_order,
                }
            )
    table = pd.DataFrame(rows).sort_values(["experiment", "target", "condition", "shRNA", "bioreplicate"])
    if table.empty:
        raise RuntimeError("No RNA-seq FASTQ pairs were found")
    return table


def build_cutrun_manifest() -> pd.DataFrame:
    require(RAW_CUTRUN_2019, "2019 CUT&RUN FASTQ root")
    rows: list[dict[str, object]] = []
    pattern19 = re.compile(r"^(?P<rep>[12])-(?P<condition>NT|shMCM3|shNONO|shPSPC1)-(?P<factor>MCM3|NONO|PSPC1)_combined_R1\.fastq\.gz$", re.I)
    for r1 in sorted(RAW_CUTRUN_2019.rglob("*_combined_R1.fastq.gz")):
        hit = pattern19.match(r1.name)
        if hit is None:
            continue
        r2 = r1.with_name(r1.name.replace("_R1.fastq.gz", "_R2.fastq.gz"))
        if not r2.is_file():
            raise FileNotFoundError(f"Missing paired read for {r1}")
        factor = hit.group("factor").upper()
        condition = hit.group("condition")
        rep = int(hit.group("rep"))
        # Each factor uses its NT pair plus the matching factor-KD pair.
        selected = condition == "NT" or condition == f"sh{factor}"
        rows.append({
            "library_id": f"2019_{condition}_rep{rep}_{factor}", "batch": "2019",
            "condition": condition, "factor": factor, "replicate": rep,
            "r1": str(r1.resolve()), "r2": str(r2.resolve()), "selected_factor_library": selected,
            "yeast_mass_pg": pd.NA,
        })

    expected_2019 = 12
    if len(rows) != expected_2019:
        raise RuntimeError(f"Expected {expected_2019} 2019 factor FASTQ pairs; found {len(rows)}")

    for batch, root in RAW_CUTRUN_2020.items():
        require(root, f"{batch} CUT&RUN FASTQ root")
        for r1 in sorted(root.rglob("*_combined_R1.fastq.gz")):
            hit = re.match(r"^(MCM3|NONO|PSPC1|IgG)-([12])_combined_R1\.fastq\.gz$", r1.name, re.I)
            if hit is None:
                continue
            r2 = r1.with_name(r1.name.replace("_R1.fastq.gz", "_R2.fastq.gz"))
            if not r2.is_file():
                raise FileNotFoundError(f"Missing paired read for {r1}")
            factor = hit.group(1)
            factor = "IgG" if factor.lower() == "igg" else factor.upper()
            if batch == "20200923" and hit.group(2) != "2":
                raise RuntimeError(f"Unexpected 20200923 sample name: {r1.name}")
            if batch == "20200929" and hit.group(2) != "1":
                raise RuntimeError(f"Unexpected 20200929 sample name: {r1.name}")
            rep = 2 if batch == "20200923" else 1
            rows.append({
                "library_id": f"2020_{factor}_{rep}", "batch": "2020", "condition": "NT",
                "factor": factor, "replicate": rep, "r1": str(r1.resolve()), "r2": str(r2.resolve()),
                "selected_factor_library": factor != "IgG", "yeast_mass_pg": 150.0,
            })
    table = pd.DataFrame(rows).sort_values(["batch", "factor", "condition", "replicate"]).reset_index(drop=True)
    expected_2020 = 8
    if int((table["batch"] == "2020").sum()) != expected_2020:
        raise RuntimeError(f"Expected {expected_2020} 2020 factor/IgG FASTQ pairs")
    return table


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()
    for directory in (RNA_ROOT / "metadata", CUTRUN_ROOT / "metadata", REFERENCE_DIR):
        directory.mkdir(parents=True, exist_ok=True)
    rna = build_rna_manifest()
    cutrun = build_cutrun_manifest()
    if not args.dry_run:
        rna.to_csv(RNA_ROOT / "metadata" / "Samples.tsv", sep="\t", index=False)
        cutrun.to_csv(CUTRUN_ROOT / "metadata" / "Samples.tsv", sep="\t", index=False)
    print(f"RNA-seq libraries: {len(rna)}; selected: {int(rna.selected.sum())}")
    print(f"CUT&RUN libraries: {len(cutrun)}; selected 2019 factor libraries: {int(((cutrun.batch == '2019') & cutrun.selected_factor_library).sum())}; 2020: {int((cutrun.batch == '2020').sum())}")


if __name__ == "__main__":
    main()
