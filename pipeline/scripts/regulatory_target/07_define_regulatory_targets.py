#!/usr/bin/env python3
"""Intersect promoter-bound genes with DEGs and render regulatory figures."""

#Before you run this script, please replace directory placeholders with your own directory

from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
from pathlib import Path

PIPELINE_ROOT = Path(__file__).resolve().parents[1]
if str(PIPELINE_ROOT) not in sys.path:
    sys.path.insert(0, str(PIPELINE_ROOT))

from config import CUTRUN_ROOT, REGULATORY_ROOT, RUN_ROOT


HELPERS = Path(__file__).resolve().parent / "figures_rendering"
VISUALS = REGULATORY_ROOT / "visuals"


def run_r(script: str) -> None:
    environment = os.environ.copy()
    environment["PATH"] = str(Path(sys.executable).parent) + os.pathsep + environment.get("PATH", "")
    executable = shutil.which("Rscript", path=environment["PATH"])
    if executable is None:
        raise FileNotFoundError("Rscript")
    for variable in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"):
        environment[variable] = "1"
    cache = REGULATORY_ROOT / "data" / ".matplotlib"
    cache.mkdir(parents=True, exist_ok=True)
    environment["MPLCONFIGDIR"] = str(cache)
    subprocess.run([str(executable), str(HELPERS / script), str(RUN_ROOT)], check=True, env=environment)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()
    if args.dry_run:
        print("[DRY-RUN] Would calculate regulatory targets and render their figures.")
        return
    promoter_input = CUTRUN_ROOT / "data" / "figure_inputs" / "promoter_gene_venn"
    if not promoter_input.is_dir():
        raise FileNotFoundError("Run CUT&RUN Step 05 before regulatory integration")
    run_r("render_direction.R")
    run_r("render_targets.R")
    print(f"[DONE] Regulatory-target figures: {VISUALS}")


if __name__ == "__main__":
    main()
