#!/usr/bin/env python3
"""Intersect promoter-bound genes with DEGs and render regulatory figures."""

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
    executable = Path(shutil.which("Rscript") or "/usr/local/bin/Rscript")
    if not executable.is_file():
        raise FileNotFoundError("Rscript")
    environment = os.environ.copy()
    environment["PATH"] = str(Path(sys.executable).parent) + os.pathsep + environment.get("PATH", "")
    for variable in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"):
        environment[variable] = "1"
    cache = REGULATORY_ROOT / "data" / ".matplotlib"
    cache.mkdir(parents=True, exist_ok=True)
    environment["MPLCONFIGDIR"] = str(cache)
    subprocess.run([str(executable), str(HELPERS / script), str(RUN_ROOT)], check=True, env=environment)


OUTPUTS = {
    "promoter_vs_DEG_direction_MCM3.png": "promoter_vs_DEG_direction_MCM3.png",
    "promoter_vs_DEG_direction_NONO.png": "promoter_vs_DEG_direction_NONO.png",
    "promoter_vs_DEG_direction_PSPC1.png": "promoter_vs_DEG_direction_PSPC1.png",
    "Fig_promoter_vs_DEG_direction.png": "promoter_vs_DEG_direction_all.png",
    "venn_regulatory_targets_MCM3_NONO_PSPC1.png": "Venn_target.png",
    "venn_regulatory_targets_MCM3_NONO_PSPC1_no_numbers.png": "Venn_target_noNumbers.png",
    "group_A_pie.png": "Pie_GroupA.png",
    "group_B_pie.png": "Pie_GroupB.png",
    "group_C_pie.png": "Pie_GroupC.png",
    "group_D_pie.png": "Pie_GroupD.png",
    "group_A_pie_no_labels.png": "Pie_GroupA_noLabels.png",
    "group_B_pie_no_labels.png": "Pie_GroupB_noLabels.png",
    "group_C_pie_no_labels.png": "Pie_GroupC_noLabels.png",
    "group_D_pie_no_labels.png": "Pie_GroupD_noLabels.png",
    "groups_A_B_C_D_pies.png": "Pie_Groups.png",
    "pie_targets_promoter_binding_NONO.png": "Pie_TargetBinding_NONO.png",
    "pie_targets_promoter_binding_PSPC1.png": "Pie_TargetBinding_PSPC1.png",
    "pie_targets_promoter_binding_NONO_no_numbers.png": "Pie_TargetBinding_NONO_noNumbers.png",
    "pie_targets_promoter_binding_PSPC1_no_numbers.png": "Pie_TargetBinding_PSPC1_noNumbers.png",
    "metaprofile_geneBody_A.single.paper.png": "Metaprofile_GroupA.png",
    "metaprofile_geneBody_B.single.paper.png": "Metaprofile_GroupB.png",
    "metaprofile_geneBody_C.single.paper.png": "Metaprofile_GroupC.png",
    "metaprofile_geneBody_D.single.paper.png": "Metaprofile_GroupD.png",
    "PanelC_metaprofile_MCM3_targets_geneBody_TSS_TES_scaled_pub.png": "Metaprofile_MCM3_target.png",
}


def publish() -> None:
    VISUALS.mkdir(parents=True, exist_ok=True)
    for source_name, target_name in OUTPUTS.items():
        source, target = VISUALS / source_name, VISUALS / target_name
        if source == target:
            if not target.is_file():
                raise FileNotFoundError(target)
            continue
        if not source.is_file():
            if target.is_file():
                continue
            raise FileNotFoundError(source)
        target.unlink(missing_ok=True)
        source.replace(target)
    for path in VISUALS.glob("._*"):
        path.unlink()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--publish-only", action="store_true")
    args = parser.parse_args()
    if args.dry_run:
        print("[DRY-RUN] Would calculate regulatory targets and render their figures.")
        return
    if not args.publish_only:
        promoter_input = CUTRUN_ROOT / "data" / "figure_inputs" / "promoter_gene_venn"
        if not promoter_input.is_dir():
            raise FileNotFoundError("Run CUT&RUN Step 05 before regulatory integration")
        run_r("render_direction.R")
        run_r("render_targets.R")
    publish()
    print(f"[DONE] Regulatory-target figures: {VISUALS}")


if __name__ == "__main__":
    main()
