#!/usr/bin/env python3
"""Step 01: prepare references, manifests, and accepted local inputs.

The raw mode builds the full reference indexes and validates original FASTQs.
The accepted mode additionally materializes the final-local accepted
intermediates used by the publication rendering workflow.
"""

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
    parser.add_argument("--threads", type=int, default=max(1, min(16, os.cpu_count() or 1)))
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()

    common = ("--threads", str(args.threads), *( ("--dry-run",) if args.dry_run else () ))
    reference_args = (*common, *( ("--full-index",) if args.source == "raw" else () ))
    run("_references.py", *reference_args)
    run("_inputs.py", *( ("--dry-run",) if args.dry_run else () ))
    if args.source == "accepted":
        run("_materialize_intermediates.py", *( ("--dry-run",) if args.dry_run else () ))
    status = "DRY-RUN" if args.dry_run else "DONE"
    print(f"[{status}] Step 01 preparation completed ({args.source} source).")


if __name__ == "__main__":
    main()
