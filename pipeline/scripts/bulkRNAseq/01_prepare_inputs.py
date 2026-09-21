#!/usr/bin/env python3
"""Prepare references, manifests, and accepted local inputs."""

#Before you run this script, please replace directory placeholders with your own directory

from __future__ import annotations

import argparse
import os
import subprocess
import sys
from pathlib import Path


HERE = Path(__file__).resolve().parent


def run(script: str, *arguments: str) -> None:
    command = [sys.executable, str(HERE / script), *arguments]
    print("[RUN]", " ".join(command), flush=True)
    subprocess.run(command, check=True)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", choices=("accepted", "raw"), default="accepted")
    parser.add_argument(
        "--threads", type=int, default=max(1, min(16, os.cpu_count() or 1))
    )
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()

    dry_run_args = ("--dry-run",) if args.dry_run else ()
    common_args = ("--threads", str(args.threads), *dry_run_args)
    reference_args = (
        *common_args,
        *(("--full-index",) if args.source == "raw" else ()),
    )
    run("_prepare_references.py", *reference_args)
    run("_build_manifests.py", *dry_run_args)
    if args.source == "accepted":
        run("_stage_accepted_inputs.py", *dry_run_args)
    status = "DRY-RUN" if args.dry_run else "DONE"
    print(f"[{status}] Step 01 preparation completed ({args.source} source).")


if __name__ == "__main__":
    main()
