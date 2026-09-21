"""Three-set CUT&RUN peak Venn renderer."""

from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib_venn import venn3, venn3_circles


COLORS = ("#A80F14", "#F39C12", "#1F4AA8")
REGION_COLORS = {
    "100": "#A80F14", "010": "#F39C12", "001": "#1F4AA8",
    "110": "#EE7D77", "101": "#8C6E93", "011": "#8A8F5A", "111": "#C9C9C9",
}
TEXT_OFFSETS = {"010": (-0.01, 0.00), "011": (0.015, -0.02)}


def plot_triple_venn(
    sets: tuple[set[str], set[str], set[str]],
    labels: tuple[str, str, str],
    output: Path,
    show_numbers: bool = True,
    show_totals: bool = True,
) -> None:
    a, b, c = sets
    values = (
        len(a - b - c), len(b - a - c), len((a & b) - c),
        len(c - a - b), len((a & c) - b), len((b & c) - a), len(a & b & c),
    )
    figure = plt.figure(figsize=(6.3, 6.3), dpi=300, facecolor="white")
    axis = figure.add_axes([0.06, 0.10, 0.88, 0.78], facecolor="white")
    diagram = venn3(subsets=values, set_labels=("", "", ""), ax=axis)
    circles = venn3_circles(subsets=values, ax=axis, linewidth=1.0, color="#222222")
    for circle in circles:
        if circle is not None:
            circle.set_alpha(0.45)
    for region, color in REGION_COLORS.items():
        patch = diagram.get_patch_by_id(region)
        if patch is not None:
            patch.set_facecolor(color)
            patch.set_alpha(0.92 if region in {"100", "010", "001"} else 0.70)
            patch.set_edgecolor("none")
            patch.set_linewidth(0.0)
            patch.set_hatch(None)
        text = diagram.get_label_by_id(region)
        if text is None:
            continue
        text.set_visible(show_numbers)
        if show_numbers:
            value = int(text.get_text().replace(",", "") or 0)
            text.set_fontsize(12)
            text.set_fontweight("bold")
            text.set_color("#F4F4F4" if region in {"100", "110", "101", "001"} and value >= 100 else "#2A2A2A")
            dx, dy = TEXT_OFFSETS.get(region, (0.0, 0.0))
            x, y = text.get_position()
            text.set_position((x + dx, y + dy))
    if show_numbers:
        nudges = {
            "100": (-0.015, 0.015), "010": (-0.030, 0.005), "001": (0.010, -0.010),
            "110": (0.000, 0.010), "101": (0.000, -0.020), "011": (-0.010, -0.006),
            "111": (0.000, -0.002),
        }
        for region, (dx, dy) in nudges.items():
            text = diagram.get_label_by_id(region)
            if text is not None:
                x, y = text.get_position()
                text.set_position((x + dx, y + dy))
    axis.set_aspect("equal", adjustable="box")
    axis.set_axis_off()
    totals = tuple(len(item) for item in sets)
    positions = ((0.10, 0.84), (0.78, 0.77), (0.76, 0.16))
    for label, total, position in zip(labels, totals, positions):
        text = f"{label}\n({total:,})" if show_totals else label
        figure.text(*position, text, ha="left", va="center", fontsize=14, fontweight="bold", color="#2A2A2A")
    output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(figure)
