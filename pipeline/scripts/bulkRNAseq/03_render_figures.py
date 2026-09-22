#!/usr/bin/env python3
"""Render all bulk RNA-seq publication figures from packaged quantifications."""

from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent
PIPELINE_ROOT = SCRIPT_DIR.parent
if str(PIPELINE_ROOT) not in sys.path:
    sys.path.insert(0, str(PIPELINE_ROOT))

from config import REFERENCE_DIR, RNA_ROOT


def executable(name: str) -> str:
    path = shutil.which(name)
    if path:
        return path
    raise FileNotFoundError(name)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()

    manifest = RNA_ROOT / "metadata" / "Salmon_quantifications.tsv"
    if not manifest.is_file() and not args.dry_run:
        raise FileNotFoundError(f"Missing packaged quantification manifest: {manifest}")

    r_command = [
        executable("Rscript"),
        str(SCRIPT_DIR / "_run_deg_analysis.R"),
        str(RNA_ROOT),
        str(REFERENCE_DIR),
    ]
    venn_command = [sys.executable, str(SCRIPT_DIR / "_render_venn.py")]
    if args.dry_run:
        print("[DRY-RUN]", " ".join(r_command))
        print("[DRY-RUN]", " ".join(venn_command))
        return

    subprocess.run(r_command, check=True)
    subprocess.run(venn_command, check=True)
    print(f"[DONE] Bulk RNA-seq figures: {RNA_ROOT / 'visuals'}")


if __name__ == "__main__":
    main()
