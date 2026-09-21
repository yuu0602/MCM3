#!/usr/bin/env python3
"""Create analysis manifests from recursively discovered raw FASTQs."""

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
    RAW_CUTRUN_DIR,
    RAW_RNA_DIR,
    REFERENCE_DIR,
    RNA_ROOT,
)


LANES = ("L001", "L002", "L003", "L004")
MCM3_SH2_LANES = frozenset({"L002", "L003"})
NONO_SELECTED = frozenset(
    {
        ("sh455", "L001"),
        ("sh455", "L003"),
        ("sh777", "L001"),
        ("sh777", "L003"),
    }
)
PSPC1_SELECTED = frozenset(
    {
        ("sh197", "L001"),
        ("sh197", "L002"),
        ("sh873", "L002"),
        ("sh873", "L003"),
    }
)
FASTQ_SUFFIXES = (".fastq.gz", ".fq.gz")
RNA_EXPERIMENTS = ("20161213", "201709", "20180718")
CUTRUN_BATCHES = ("20200923", "20200929")


def require(path: Path, label: str) -> None:
    if not path.exists():
        raise FileNotFoundError(f"Missing {label}: {path}")


def normalized_stem(path: Path) -> str:
    """Return a filename-only stem for supported compressed FASTQ names."""
    for suffix in FASTQ_SUFFIXES:
        if path.name.endswith(suffix):
            return path.name.removesuffix(suffix).replace("-", "_")
    raise ValueError(f"Unsupported FASTQ suffix: {path}")


def paired_fastq(r1: Path, read1_marker: str, read2_marker: str) -> Path:
    """Locate the mate beside a recursively discovered read-one FASTQ."""
    for suffix in FASTQ_SUFFIXES:
        read1_suffix = f"{read1_marker}{suffix}"
        if r1.name.endswith(read1_suffix):
            return r1.with_name(
                r1.name.removesuffix(read1_suffix) + f"{read2_marker}{suffix}"
            )
    raise ValueError(f"Read-one filename does not end in {read1_marker}: {r1}")


def recursive_read_ones(root: Path, marker: str) -> list[Path]:
    """Find read-one FASTQs recursively using filename patterns only."""
    return sorted(
        {
            path
            for suffix in FASTQ_SUFFIXES
            for path in root.rglob(f"*{marker}{suffix}")
            if path.is_file()
        }
    )


def metadata_token(path: Path, choices: tuple[str, ...]) -> str | None:
    """Read a known cohort token from the filename, then its parent directories."""
    pattern = re.compile(
        rf"(?:^|[_-])({'|'.join(map(re.escape, choices))})(?:[_-]|$)"
    )
    for value in (path.name, *(parent.name for parent in path.parents)):
        hit = pattern.search(value)
        if hit:
            return hit.group(1)
    return None


def lane_from_name(name: str) -> str:
    hit = re.search(r"_L(\d{3})_", name)
    if hit is None:
        raise ValueError(f"Cannot infer L00X biological replicate from {name}")
    return f"L{hit.group(1)}"


def classify_rna(name: str) -> tuple[str, str, str, str]:
    text = normalized_stem(Path(name)).lower()
    hit = re.search(r"(?:^|_)(20161213|201709|20180718)(?:_|$)", text)
    if hit is None:
        raise ValueError(f"RNA-seq filename lacks an experiment date: {name}")
    experiment = hit.group(1)
    if experiment == "20161213":
        if "mcm3kd_nt" in text:
            return experiment, "NT", "MCM3", "NT"
        if re.search(r"(?:^|_)sh1(?:_|$)", text):
            return experiment, "KD", "MCM3", "sh1"
        if re.search(r"(?:^|_)sh2(?:_|$)", text):
            return experiment, "KD", "MCM3", "sh2"
    if experiment == "201709":
        if re.search(r"(?:^|_)nt(?:_|$)", text):
            return experiment, "NT", "MCM3", "NT"
        if "sh3_6" in text:
            return experiment, "KD", "MCM3", "sh3_6"
        if "sh3_10" in text:
            return experiment, "KD", "MCM3", "sh3_10"
    if experiment == "20180718":
        if re.search(r"(?:^|_)nt(?:_|$)", text):
            return experiment, "NT", "NT", "NT"
        hit = re.search(r"(?:^|_)(nono|pspc1)_sh(\d+)(?:_|$)", text)
        if hit:
            return experiment, "KD", hit.group(1).upper(), f"sh{hit.group(2)}"
    raise ValueError(f"Unrecognized RNA-seq sample name: {name}")


def include_rna(
    experiment: str,
    condition: str,
    target: str,
    shrna: str,
    lane: str,
) -> tuple[bool, str]:
    """Encode the final RNA-seq sample selection."""
    if target == "MCM3":
        keep = experiment == "20161213" and (
            (condition == "NT" and lane in LANES)
            or (shrna == "sh1" and lane in LANES)
            or (shrna == "sh2" and lane in MCM3_SH2_LANES)
        )
        return keep, "MCM3_NT_vs_sh1_sh2" if keep else "not_selected"
    if target in {"NONO", "PSPC1"}:
        if experiment != "20180718":
            return False, "not_selected"
        if condition == "NT":
            keep = lane in LANES
        elif target == "NONO":
            keep = (shrna, lane) in NONO_SELECTED
        else:
            keep = (shrna, lane) in PSPC1_SELECTED
        return keep, f"{target}_NT_vs_pooledKD" if keep else "not_selected"
    if target == "NT" and experiment == "20180718":
        return lane in LANES, "shared_20180718_NT"
    return False, "not_selected"


def legacy_nono_pspc1_context(
    experiment: str, condition: str, target: str, shrna: str, lane: str
) -> tuple[bool, int | None]:
    """Return the fixed NONO/PSPC1 normalization cohort and its order.

    NONO and PSPC1 were accepted from the original per-target analysis.  That
    analysis calculated TMM factors and its expression filter on this fixed
    22-library cohort before fitting the two 20180718 target contrasts.  The
    MCM3 contrast itself is still fitted independently from the current
    sh1+sh2 sample set.
    """
    key = (experiment, condition, target, shrna, lane)
    order: list[tuple[str, str, str, str, str]] = []
    order.extend(("20161213", "KD", "MCM3", "sh1", value) for value in LANES)
    order.extend(("201709", "KD", "MCM3", "sh3_6", value) for value in ("L001", "L003"))
    order.extend(("20161213", "NT", "MCM3", "NT", value) for value in ("L002", "L004"))
    order.extend(("201709", "NT", "MCM3", "NT", value) for value in ("L001", "L002"))
    order.extend(
        ("20180718", "KD", "NONO", "sh455", value)
        for value in ("L001", "L003")
    )
    order.extend(
        ("20180718", "KD", "NONO", "sh777", value)
        for value in ("L001", "L003")
    )
    order.extend(("20180718", "NT", "NT", "NT", value) for value in LANES)
    order.extend(
        ("20180718", "KD", "PSPC1", "sh197", value)
        for value in ("L001", "L002")
    )
    order.extend(
        ("20180718", "KD", "PSPC1", "sh873", value)
        for value in ("L002", "L003")
    )
    try:
        return True, order.index(key) + 1
    except ValueError:
        return False, None


def build_rna_manifest() -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    require(RAW_RNA_DIR, "RNA-seq FASTQ directory")
    for r1 in recursive_read_ones(RAW_RNA_DIR, "_R1_001"):
        if r1.name.startswith(("._", "Undetermined", "Unindexed_Reads")):
            continue
        r2 = paired_fastq(r1, "_R1_001", "_R2_001")
        if not r2.is_file():
            raise FileNotFoundError(f"Missing paired read for {r1}")
        experiment = metadata_token(r1, RNA_EXPERIMENTS)
        if experiment is None:
            raise ValueError(f"RNA-seq sample lacks an experiment date: {r1}")
        # Classification remains filename-based; the parent fallback supplies
        # only the legacy experiment identifier when nested layouts omit it.
        experiment, condition, target, shrna = classify_rna(
            f"{experiment}_{r1.name}"
        )
        lane = lane_from_name(r1.name)
        selected, analysis = include_rna(experiment, condition, target, shrna, lane)
        legacy_context, legacy_order = legacy_nono_pspc1_context(
            experiment, condition, target, shrna, lane
        )
        rows.append(
            {
                "sample_id": f"{experiment}_{target}_{condition}_{shrna}_{lane}",
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
    table = pd.DataFrame(rows).sort_values(
        ["experiment", "target", "condition", "shRNA", "bioreplicate"]
    )
    if table.empty:
        raise RuntimeError("No RNA-seq FASTQ pairs were found")
    return table


def build_cutrun_manifest() -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    require(RAW_CUTRUN_DIR, "CUT&RUN FASTQ directory")
    for r1 in recursive_read_ones(RAW_CUTRUN_DIR, "_R1"):
        name = normalized_stem(r1)
        batch_id = metadata_token(r1, CUTRUN_BATCHES)
        factor = re.search(r"(?:^|_)(MCM3|NONO|PSPC1|IgG)(?:_|$)", name, re.I)
        replicate = re.search(r"(?:^|_)([12])_R1$", name)
        if batch_id is None or factor is None or replicate is None:
            continue
        r2 = paired_fastq(r1, "_R1", "_R2")
        if not r2.is_file():
            raise FileNotFoundError(f"Missing paired read for {r1}")
        factor_id = factor.group(1)
        factor_id = "IgG" if factor_id.lower() == "igg" else factor_id.upper()
        rep = int(replicate.group(1))
        if (batch_id, rep) not in {("20200923", 2), ("20200929", 1)}:
            raise RuntimeError(f"Unexpected CUT&RUN sample name: {r1.name}")
        rows.append(
            {
                "library_id": f"{factor_id}-{rep}",
                "batch": "2020",
                "condition": "NT",
                "factor": factor_id,
                "replicate": rep,
                "r1": str(r1.resolve()),
                "r2": str(r2.resolve()),
                "selected_factor_library": factor_id != "IgG",
                "yeast_mass_pg": 150.0,
            }
        )
    table = pd.DataFrame(rows).sort_values(
        ["batch", "factor", "condition", "replicate"]
    ).reset_index(drop=True)
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
    factor_libraries = int(cutrun.selected_factor_library.sum())
    igg_libraries = int((cutrun.factor == "IgG").sum())
    print(
        f"CUT&RUN libraries: {len(cutrun)}; "
        f"2020 factor libraries: {factor_libraries}; "
        f"matched IgG libraries: {igg_libraries}"
    )


if __name__ == "__main__":
    main()
