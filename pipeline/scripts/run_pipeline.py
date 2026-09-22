#!/usr/bin/env python3
"""Run the MCM3 publication workflow in dependency order."""

from __future__ import annotations

import argparse
import os
import subprocess
import sys
from pathlib import Path


PREPARE = "bulkRNAseq/01_prepare_inputs.py"
RNASEQ = "bulkRNAseq/02_quantify_and_test.py"
CUTRUN = "CUT&RUN/04_align_and_normalize.py"
BINDING = "CUT&RUN/05_call_peaks_and_define_binding.py"
CUTRUN_FIGURES = "CUT&RUN/06_render_figures.py"
REGULATORY = "regulatory_target/07_define_regulatory_targets.py"


def workflow(source: str) -> dict[str, tuple[str, ...]]:
    raw_cutrun = (CUTRUN,) if source == "raw" else ()
    figures = (PREPARE, RNASEQ, *raw_cutrun, BINDING, CUTRUN_FIGURES, REGULATORY)
    return {
        "prepare": (PREPARE,),
        "rna": (PREPARE, RNASEQ),
        "cutrun": (PREPARE, *raw_cutrun, BINDING, CUTRUN_FIGURES),
        "regulatory": figures,
        "figures": figures,
        "all": figures,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--step", choices=("all", "prepare", "rna", "cutrun", "regulatory", "figures"), default="all")
    parser.add_argument("--source", choices=("accepted", "raw"), default="accepted")
    parser.add_argument("--threads", type=int, default=max(1, min(16, os.cpu_count() or 1)))
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()

    root = Path(__file__).resolve().parent
    python = os.environ.get("MCM3_FINAL_PYTHON") or sys.executable

    for relative in workflow(args.source)[args.step]:
        command = [python, str(root / relative)]
        if relative in {PREPARE, RNASEQ, CUTRUN}:
            command += ["--threads", str(args.threads)]
        if relative == PREPARE:
            command += ["--source", args.source]
        if args.source == "raw" and relative == RNASEQ:
            command.append("--requantify")
        if args.source == "raw" and relative in {CUTRUN, BINDING}:
            command.append("--from-raw")
        if args.dry_run:
            command.append("--dry-run")
        print("[RUN]", " ".join(command), flush=True)
        subprocess.run(command, check=True)


if __name__ == "__main__":
    main()
