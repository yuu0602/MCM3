#!/usr/bin/env python3
"""Generate RNA-seq DEG outputs from local Salmon quantifications.

The default path uses the accepted Salmon products materialized by
``01_prepare.py``.  ``--requantify`` is the full raw-FASTQ
path and rebuilds every selected quantification before the same DEG script.
"""

from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pandas as pd

PIPELINE_ROOT = Path(__file__).resolve().parents[1]
if str(PIPELINE_ROOT) not in sys.path:
    sys.path.insert(0, str(PIPELINE_ROOT))

from config import GENCODE_TRANSCRIPTS_NAME, REFERENCE_DIR, RNA_ROOT


def tool(name: str, candidates: tuple[str, ...]) -> str:
    for candidate in candidates:
        resolved = shutil.which(candidate)
        if resolved:
            return resolved
        if Path(candidate).is_file():
            return candidate
    raise FileNotFoundError(f"Required executable not found: {name}")


def render_venn(dry_run: bool) -> None:
    """Render the final RNA-seq Venn panels as part of Step 02."""
    command = [sys.executable, str(Path(__file__).with_name("_render_rnaseq_venn.py"))]
    if dry_run:
        command.append("--dry-run")
    print("[RUN]", " ".join(command))
    subprocess.run(command, check=True)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--threads", type=int, default=max(1, min(16, os.cpu_count() or 1)))
    parser.add_argument("--requantify", action="store_true", help="Rebuild Salmon quantifications from original FASTQs")
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()

    manifest_path = RNA_ROOT / "metadata" / "Samples.tsv"
    if not manifest_path.is_file():
        raise FileNotFoundError("Run 01_prepare.py before RNA-seq quantification")
    rscript = tool("Rscript", ("Rscript", "/usr/local/bin/Rscript"))
    quant_manifest = RNA_ROOT / "metadata" / "Salmon_quantifications.tsv"
    if not args.requantify:
        if not quant_manifest.is_file() and not args.dry_run:
            raise FileNotFoundError(
                "Missing materialized Salmon_quantifications.tsv. Run 01_prepare.py "
                "or use --requantify for the raw-FASTQ path."
            )
        command = [rscript, str(Path(__file__).with_name("_rnaseq_deg.R")), str(RNA_ROOT), str(REFERENCE_DIR)]
        print("[RUN]", " ".join(command))
        if not args.dry_run:
            subprocess.run(command, check=True)
        render_venn(args.dry_run)
        return

    transcriptome = REFERENCE_DIR / GENCODE_TRANSCRIPTS_NAME
    index = REFERENCE_DIR / "salmon_index_M25"
    salmon = tool("salmon", ("salmon", "/opt/homebrew/bin/salmon", "/opt/anaconda3/bin/salmon"))
    manifest = pd.read_csv(manifest_path, sep="\t", dtype=str)
    selected_mask = manifest["selected"].str.lower().eq("true")
    legacy_mask = manifest.get(
        "legacy_nono_pspc1_context", pd.Series(False, index=manifest.index)
    ).astype(str).str.lower().eq("true")
    selected = manifest.loc[selected_mask | legacy_mask].copy()
    if selected.empty:
        raise RuntimeError("RNA-seq manifest has no selected samples")
    if not args.dry_run and not transcriptome.is_file():
        raise FileNotFoundError(
            f"Missing {transcriptome}. Run 01_prepare.py before RNA-seq quantification."
        )

    if not index.is_dir():
        command = [salmon, "index", "-t", str(transcriptome), "-i", str(index), "-k", "31", "-p", str(args.threads)]
        print("[RUN]", " ".join(command))
        if not args.dry_run:
            subprocess.run(command, check=True)

    quant_root = RNA_ROOT / "salmon"
    rows: list[dict[str, str]] = []
    for record in selected.to_dict("records"):
        out_dir = quant_root / str(record["sample_id"])
        quant = out_dir / "quant.sf"
        command = [
            salmon, "quant", "-i", str(index), "-l", "A", "-1", str(record["r1"]), "-2", str(record["r2"]),
            "-p", str(args.threads), "--validateMappings", "-o", str(out_dir),
        ]
        print("[RUN]", " ".join(command))
        if not args.dry_run and not quant.is_file():
            subprocess.run(command, check=True)
        rows.append({**record, "quant_sf": str(quant)})

    if not args.dry_run:
        pd.DataFrame(rows).to_csv(quant_manifest, sep="\t", index=False)
        subprocess.run([rscript, str(Path(__file__).with_name("_rnaseq_deg.R")), str(RNA_ROOT), str(REFERENCE_DIR)], check=True)
    render_venn(args.dry_run)


if __name__ == "__main__":
    main()
