#!/usr/bin/env python3
"""Copy accepted computational intermediates into the final run directory.

This is the default, time-efficient execution path.  It makes the final run
self-contained before any result table or publication figure is rendered.  It
does not copy a prior figure or DEG result: later stages generate those from
the local Salmon quantifications, accepted peak evidence, and bigWig tracks.
"""

from __future__ import annotations

import argparse
import os
import shutil
import sys
from pathlib import Path

import pandas as pd

PIPELINE_ROOT = Path(__file__).resolve().parents[1]
if str(PIPELINE_ROOT) not in sys.path:
    sys.path.insert(0, str(PIPELINE_ROOT))

from config import ACCEPTED_INTERMEDIATE_SOURCES, CUTRUN_ROOT, GENCODE_GTF_NAME, REFERENCE_DIR, RNA_ROOT, RUN_ROOT


DESTINATIONS = {
    "reference_gtf": REFERENCE_DIR / GENCODE_GTF_NAME,
    "rna_salmon": RNA_ROOT / "salmon",
    "cutrun_2019_bowtie2": CUTRUN_ROOT / "01_bowtie2" / "2019",
    "cutrun_2019_bam": CUTRUN_ROOT / "02_bam" / "2019",
    "cutrun_2019_peaks": CUTRUN_ROOT / "04_peaks" / "2019",
    "cutrun_2020_bowtie2": CUTRUN_ROOT / "01_bowtie2" / "2020",
    "cutrun_2020_bam": CUTRUN_ROOT / "02_bam" / "2020",
    "cutrun_2020_peaks": CUTRUN_ROOT / "04_peaks" / "2020",
    "cutrun_bigwig": CUTRUN_ROOT / "03_bigwig",
    "cutrun_membership": CUTRUN_ROOT / "data" / "figure_inputs",
}

IGV_TRACKS = {
    "MCM3_six_library_batchBalanced_mean.bw": "MCM3_mean.bw",
    "NONO_six_library_batchBalanced_mean.bw": "NONO_mean.bw",
    "PSPC1_six_library_batchBalanced_mean.bw": "PSPC1_mean.bw",
    "IgG_2020_to_2019_commonBridge_mean.bw": "IgG_mean.bw",
}

BIGWIG_ROOT_FILES = {
    "batch_adjustment_parameters.tsv": "IGV_representation_batch_adjustment.tsv",
    "six_library_track_manifest.tsv": "IGV_representation_libraries.tsv",
}

SALMON_ROOT_FILES = {
    "samples.tsv": "Samples_all.tsv",
    "samples_with_quant.tsv": "Quantifications_all.tsv",
    "biorep_samples_with_quant.tsv": "Quantifications_bioreplicates.tsv",
    "biorep_qc_flags.tsv": "QC_bioreplicates_flags.tsv",
    "biorep_qc_summary.tsv": "QC_bioreplicates_summary.tsv",
    "lane_biorep_qc_flags.tsv": "QC_lane_bioreplicates_flags.tsv",
    "lane_biorep_qc_summary.tsv": "QC_lane_bioreplicates_summary.tsv",
    "lane_qc_flags.tsv": "QC_lanes_flags.tsv",
    "lane_qc_summary.tsv": "QC_lanes_summary.tsv",
    "replicate_qc_flags.tsv": "QC_replicates_flags.tsv",
    "replicate_qc_summary.tsv": "QC_replicates_summary.tsv",
}


def ignored(path: Path) -> bool:
    return path.name == ".DS_Store" or path.name.startswith("._") or path.name.endswith(".bak") or "__pycache__" in path.parts


def copy_tree(source: Path, destination: Path, dry_run: bool) -> tuple[int, int]:
    """Copy only missing or size-changed files and return file/byte counts."""
    if source.is_file():
        matches = destination.is_file() and destination.stat().st_size == source.stat().st_size
        if matches:
            return 0, 0
        print(f"[COPY] {source} -> {destination}", flush=True)
        if not dry_run:
            destination.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(source, destination)
        return 1, source.stat().st_size
    if not source.is_dir():
        raise FileNotFoundError(f"Missing accepted intermediate source: {source}")
    copied = 0
    copied_bytes = 0
    for root, dirs, files in os.walk(source):
        root_path = Path(root)
        dirs[:] = [name for name in dirs if not ignored(root_path / name)]
        relative = root_path.relative_to(source)
        target_dir = destination / relative
        for name in files:
            current = root_path / name
            if ignored(current):
                continue
            target = target_dir / name
            size_matches = target.is_file() and target.stat().st_size == current.stat().st_size
            if size_matches:
                continue
            copied += 1
            copied_bytes += current.stat().st_size
            print(f"[COPY] {current} -> {target}", flush=True)
            if not dry_run:
                target.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(current, target)
    return copied, copied_bytes


def copy_bigwig_tree(source: Path, destination: Path, dry_run: bool) -> tuple[int, int]:
    """Materialize accepted bigWigs in the final IGV/intermediate layout."""
    if not source.is_dir():
        raise FileNotFoundError(f"Missing accepted intermediate source: {source}")
    copied = 0
    copied_bytes = 0
    for current in source.rglob("*"):
        if not current.is_file() or ignored(current):
            continue
        relative = current.relative_to(source)
        if len(relative.parts) == 1 and current.name in BIGWIG_ROOT_FILES:
            target = destination / BIGWIG_ROOT_FILES[current.name]
        elif relative.parts[0] == "replicate_mean":
            if len(relative.parts) >= 2 and relative.parts[1] == "display10x":
                if current.name not in IGV_TRACKS:
                    continue
                target = destination / "IGV_representation" / IGV_TRACKS[current.name]
            else:
                target = destination / "mean_intermediates" / Path(*relative.parts[1:])
        else:
            target = destination / relative
        if target.is_file() and target.stat().st_size == current.stat().st_size:
            continue
        copied += 1
        copied_bytes += current.stat().st_size
        print(f"[COPY] {current} -> {target}", flush=True)
        if not dry_run:
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(current, target)

    source_manifest = source / "replicate_mean" / "display10x" / "display10x_manifest.tsv"
    if source_manifest.is_file() and not dry_run:
        manifest = pd.read_csv(source_manifest, sep="\t")
        manifest["source_track"] = manifest["source_track"].map(
            lambda name: f"mean_intermediates/{Path(name).name}"
        )
        manifest["igv_track"] = manifest.pop("display_track").map(
            lambda name: IGV_TRACKS.get(Path(name).name, Path(name).name)
        )
        columns = ["source_track", "igv_track", *[column for column in manifest if column not in {"source_track", "igv_track"}]]
        output = destination / "IGV_representation" / "IGV_tracks_manifest.tsv"
        output.parent.mkdir(parents=True, exist_ok=True)
        manifest.loc[:, columns].to_csv(output, sep="\t", index=False)
    return copied, copied_bytes


def copy_salmon_tree(source: Path, destination: Path, dry_run: bool) -> tuple[int, int]:
    """Copy native Salmon outputs and normalize project-level table names."""
    if not source.is_dir():
        raise FileNotFoundError(f"Missing accepted intermediate source: {source}")
    copied = 0
    copied_bytes = 0
    for current in source.rglob("*"):
        if not current.is_file() or ignored(current):
            continue
        relative = current.relative_to(source)
        target_name = SALMON_ROOT_FILES.get(current.name, current.name) if len(relative.parts) == 1 else current.name
        target = destination / relative.parent / target_name
        if target.is_file() and target.stat().st_size == current.stat().st_size:
            continue
        copied += 1
        copied_bytes += current.stat().st_size
        print(f"[COPY] {current} -> {target}", flush=True)
        if not dry_run:
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(current, target)
    return copied, copied_bytes


def copy_cutrun_membership(source: Path, destination: Path, dry_run: bool) -> tuple[int, int]:
    """Copy only the accepted peak/promoter inputs required by downstream stages."""
    tag = "p1e4_q5e2_fe3_union4of6"
    mappings: list[tuple[Path, Path]] = []
    for factor in ("MCM3", "NONO", "PSPC1"):
        mappings.extend((
            (
                source / "promoter_gene_venn" / f"FigX_promoter_genes_{factor}__{tag}.tsv",
                destination / "promoter_gene_venn" / f"Venn_PromoterGenes_{factor}.tsv",
            ),
            (
                source / "promoter_gene_venn" / f"FigX_promoter_gene_library_support_{factor}__{tag}.tsv",
                destination / "promoter_gene_venn" / f"Venn_PromoterGenes_{factor}_support.tsv",
            ),
            (
                source / "batch_stratified_peaks" / f"{factor}_pooled_six_library_support4of6.bed",
                destination / "batch_stratified_peaks" / f"Venn_Peaks_{factor}_support4of6.bed",
            ),
        ))
    mappings.extend((
        (
            source / "promoter_gene_venn" / f"FigX_promoter_gene_overlap_regions__{tag}.tsv",
            destination / "promoter_gene_venn" / "Venn_PromoterGenes_regions.tsv",
        ),
        (
            source / "promoter_gene_venn" / "analysis_parameters.tsv",
            destination / "promoter_gene_venn" / "AnalysisParameters.tsv",
        ),
        (
            source / "batch_stratified_peak_summary.tsv",
            destination / "Venn_Peaks_support.tsv",
        ),
        (
            source / "peak_distribution" / "reanalysis20_layout_intermediate" / "IgG_2020_union_strict_peaks.bed",
            destination / "peak_distribution" / "Peaks_IgG.bed",
        ),
    ))
    copied = 0
    copied_bytes = 0
    for current, target in mappings:
        count, size = copy_tree(current, target, dry_run)
        copied += count
        copied_bytes += size
    return copied, copied_bytes


def salmon_relative(path: str) -> Path:
    """Obtain the path below ``deg_work/salmon`` from historic metadata."""
    parts = Path(path).parts
    for index in range(len(parts) - 2):
        if parts[index:index + 2] == ("deg_work", "salmon"):
            return Path(*parts[index + 2:])
    raise ValueError(f"Could not locate deg_work/salmon within source path: {path}")


def make_local_quant_manifest(dry_run: bool) -> int:
    selected_path = RNA_ROOT / "metadata" / "Samples.tsv"
    salmon_meta = RNA_ROOT / "salmon" / "Quantifications_all.tsv"
    if not selected_path.is_file():
        raise FileNotFoundError("Run 01_prepare.py before materializing Salmon quantifications")
    if not salmon_meta.is_file():
        raise FileNotFoundError(f"Missing copied Salmon sample table: {salmon_meta}")
    selected = pd.read_csv(selected_path, sep="\t", dtype=str)
    selected_mask = selected["selected"].str.lower().eq("true")
    legacy_mask = selected.get(
        "legacy_nono_pspc1_context", pd.Series(False, index=selected.index)
    ).astype(str).str.lower().eq("true")
    selected = selected.loc[selected_mask | legacy_mask].copy()
    source = pd.read_csv(salmon_meta, sep="\t", dtype=str)
    source["r1_basename"] = source["r1_path"].map(lambda value: Path(value).name)
    selected["r1_basename"] = selected["r1"].map(lambda value: Path(value).name)
    duplicate = source.loc[source.r1_basename.duplicated(keep=False), "r1_basename"].unique()
    if len(duplicate):
        raise RuntimeError(f"Ambiguous accepted Salmon entries: {', '.join(duplicate[:5])}")
    merged = selected.merge(
        source[["r1_basename", "replicate_id", "quant_sf"]], on="r1_basename", how="left", validate="one_to_one"
    )
    if merged["quant_sf"].isna().any():
        missing = merged.loc[merged.quant_sf.isna(), "r1_basename"].tolist()
        raise RuntimeError(f"Selected FASTQs lack accepted Salmon quantification: {missing}")
    merged["quant_sf"] = merged["quant_sf"].map(lambda value: str(RNA_ROOT / "salmon" / salmon_relative(value)))
    if not all(Path(path).is_file() for path in merged["quant_sf"]):
        missing = [path for path in merged["quant_sf"] if not Path(path).is_file()]
        raise FileNotFoundError(f"Materialized quant.sf file(s) missing: {missing[:3]}")
    out = RNA_ROOT / "metadata" / "Salmon_quantifications.tsv"
    columns = [
        "sample_id", "experiment", "condition", "target", "shRNA", "bioreplicate", "r1", "r2",
        "selected", "analysis", "legacy_nono_pspc1_context", "legacy_context_order",
        "replicate_id", "quant_sf",
    ]
    if dry_run:
        print(f"[DRY-RUN] Would write local quantification manifest for {len(merged)} selected RNA-seq samples: {out}")
    else:
        merged.loc[:, columns].to_csv(out, sep="\t", index=False)
    return len(merged)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()
    rows = []
    for label, source in ACCEPTED_INTERMEDIATE_SOURCES.items():
        destination = DESTINATIONS[label]
        if args.dry_run:
            print(f"[DRY-RUN] Would materialize {label}: {source} -> {destination}")
            count = size = 0
        else:
            if label == "cutrun_bigwig":
                count, size = copy_bigwig_tree(source, destination, False)
            elif label == "rna_salmon":
                count, size = copy_salmon_tree(source, destination, False)
            elif label == "cutrun_membership":
                count, size = copy_cutrun_membership(source, destination, False)
            else:
                count, size = copy_tree(source, destination, False)
        rows.append({
            "intermediate": label,
            "accepted_source": str(source),
            "local_copy": str(destination.relative_to(RUN_ROOT)),
            "files_copied_this_run": count,
            "bytes_copied_this_run": size,
        })
    selected = make_local_quant_manifest(args.dry_run) if not args.dry_run else 26
    if not args.dry_run:
        output = RUN_ROOT / "docs" / "accepted_intermediate_manifest.tsv"
        pd.DataFrame(rows).to_csv(output, sep="\t", index=False)
    if args.dry_run:
        print(f"[DRY-RUN] Local RNA-seq quantification manifest would contain {selected} samples.")
    else:
        print(f"[DONE] Accepted intermediates materialized locally; RNA-seq quantification manifest has {selected} samples.")


if __name__ == "__main__":
    main()
