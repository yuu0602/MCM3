#!/usr/bin/env python3
"""Render the final promoter-gene Venn design."""

from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib_venn import venn3, venn3_circles


FACTORS = ("MCM3", "NONO", "PSPC1")
PLOT_TITLE = "CUT&RUN promoter-gene overlap"
REGION_COLORS = {
    "100": "#A80F14",
    "010": "#F39C12",
    "001": "#1F4AA8",
    "110": "#EE7D77",
    "101": "#8C6E93",
    "011": "#8A8F5A",
    "111": "#C9C9C9",
}


def venn_regions(a: set[str], b: set[str], c: set[str]) -> dict[str, set[str]]:
    return {
        "100_only_MCM3": a - b - c,
        "010_only_NONO": b - a - c,
        "001_only_PSPC1": c - a - b,
        "110_MCM3_NONO": (a & b) - c,
        "101_MCM3_PSPC1": (a & c) - b,
        "011_NONO_PSPC1": (b & c) - a,
        "111_MCM3_NONO_PSPC1": a & b & c,
    }


def plot_venn(sets: dict[str, set[str]], out_png: Path, show_numbers: bool) -> None:
    a, b, c = (sets[factor] for factor in FACTORS)
    regions = venn_regions(a, b, c)
    values = tuple(
        len(regions[name])
        for name in (
            "100_only_MCM3", "010_only_NONO", "110_MCM3_NONO",
            "001_only_PSPC1", "101_MCM3_PSPC1", "011_NONO_PSPC1",
            "111_MCM3_NONO_PSPC1",
        )
    )

    figure = plt.figure(figsize=(6.3, 6.3), dpi=300, facecolor="white")
    axes = figure.add_axes([0.06, 0.10, 0.88, 0.78], facecolor="white")
    venn = venn3(subsets=values, set_labels=("", "", ""), ax=axes)
    for circle in venn3_circles(subsets=values, ax=axes, linewidth=1.0, color="#222222"):
        if circle is not None:
            circle.set_alpha(0.45)
    for region_id, color in REGION_COLORS.items():
        patch = venn.get_patch_by_id(region_id)
        if patch is not None:
            patch.set_facecolor(color)
            patch.set_alpha(0.92 if region_id in {"100", "010", "001"} else 0.70)
            patch.set_edgecolor("none")
        label = venn.get_label_by_id(region_id)
        if label is not None:
            label.set_visible(show_numbers)
            if show_numbers:
                label.set_fontsize(11)
                label.set_fontweight("bold")
                label.set_color("#F4F4F4" if region_id in {"100", "110", "101", "001"} else "#2A2A2A")

    axes.set_aspect("equal", adjustable="box")
    axes.set_axis_off()
    figure.text(0.5, 0.95, PLOT_TITLE, ha="center", va="center", fontsize=15, fontweight="bold", color="#222222")
    figure.text(0.10, 0.84, f"MCM3\n({len(a):,})", ha="left", va="center", fontsize=14, fontweight="bold")
    figure.text(0.78, 0.77, f"NONO\n({len(b):,})", ha="left", va="center", fontsize=14, fontweight="bold")
    figure.text(0.76, 0.16, f"PSPC1\n({len(c):,})", ha="left", va="center", fontsize=14, fontweight="bold")
    figure.savefig(out_png, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(figure)
