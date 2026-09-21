#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt


def die(msg: str) -> None:
    raise SystemExit(msg)


def _hex_to_rgb01(h: str):
    h = h.lstrip("#")
    return (int(h[0:2], 16) / 255.0, int(h[2:4], 16) / 255.0, int(h[4:6], 16) / 255.0)


def _rgb01_to_hex(rgb):
    r, g, b = rgb
    return "#{:02x}{:02x}{:02x}".format(int(round(r * 255)), int(round(g * 255)), int(round(b * 255)))


def blend_hex(c1: str, c2: str, w: float = 0.5) -> str:
    a = _hex_to_rgb01(c1)
    b = _hex_to_rgb01(c2)
    return _rgb01_to_hex((a[0] * (1 - w) + b[0] * w, a[1] * (1 - w) + b[1] * w, a[2] * (1 - w) + b[2] * w))


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--out", required=True)
    p.add_argument("--title", required=True)
    p.add_argument("--subtitle", default="")
    p.add_argument("--a-name", required=True)
    p.add_argument("--b-name", required=True)
    p.add_argument("--c-name", required=True)
    p.add_argument("--a-total", required=True, type=int)
    p.add_argument("--b-total", required=True, type=int)
    p.add_argument("--c-total", required=True, type=int)
    for rid in ("100", "010", "001", "110", "101", "011", "111"):
        p.add_argument(f"--n{rid}", required=True, type=int)
    p.add_argument("--hide-numbers", action="store_true")
    return p.parse_args()


def main() -> None:
    args = parse_args()
    try:
        from matplotlib_venn import venn3, venn3_circles
    except Exception as e:
        die(f"Missing matplotlib-venn: {e}")

    subsets = (
        args.n100,  # A only
        args.n010,  # B only
        args.n110,  # AB only
        args.n001,  # C only
        args.n101,  # AC only
        args.n011,  # BC only
        args.n111,  # ABC
    )

    fig = plt.figure(figsize=(6.6, 6.6), dpi=300, facecolor="#FFFFFF")
    ax = fig.add_axes([0.05, 0.09, 0.90, 0.80], facecolor="#FFFFFF")

    v = venn3(subsets=subsets, set_labels=("", "", ""), ax=ax)
    circles = venn3_circles(subsets=subsets, ax=ax, linewidth=1.1, color="#555555")
    for c in circles:
        if c is not None:
            c.set_alpha(0.65)

    colA = "#A80F14"  # MCM3
    colB = "#F39C12"  # NONO
    colC = "#1F4AA8"  # PSPC1
    region_colors = {
        "100": colA,
        "010": colB,
        "001": colC,
        "110": "#EE7D77",  # A∩B
        "101": "#8C6E93",  # A∩C
        "011": "#8A8F5A",  # B∩C (olive)
        "111": "#C9C9C9",  # A∩B∩C
    }
    region_alpha = {"100": 0.92, "010": 0.82, "001": 0.88, "110": 0.72, "101": 0.72, "011": 0.72, "111": 0.95}

    for rid in ("100", "010", "001", "110", "101", "011", "111"):
        patch = v.get_patch_by_id(rid)
        if patch is None:
            continue
        patch.set_facecolor(region_colors[rid])
        patch.set_alpha(region_alpha[rid])
        patch.set_edgecolor("none")
        patch.set_linewidth(0)

    # Make all counts bold, consistent, and readable.
    label_color = {
        "100": "#F4F4F4",
        "010": "#F4F4F4",
        "001": "#F4F4F4",
        "110": "#111111",
        "101": "#F4F4F4",
        "011": "#111111",
        "111": "#111111",
    }
    nudges = {
        # Singles
        "100": (0.020, -0.020),   # MCM3 only: down-right
        "010": (-0.035, -0.020),  # NONO only: slightly more left
        "001": (-0.035, 0.020),   # PSPC1 only: slightly more left
        # Pair overlaps
        "110": (0.018, 0.010),    # MCM3∩NONO: right
        "101": (0.020, -0.008),   # MCM3∩PSPC1: right
        "011": (-0.018, 0.012),   # NONO∩PSPC1: up-left
        "111": (0.010, -0.002),
    }
    for rid in ("100", "010", "001", "110", "101", "011", "111"):
        txt = v.get_label_by_id(rid)
        if txt is None:
            continue
        if args.hide_numbers:
            txt.set_visible(False)
            continue
        if rid in {"100", "010", "001"}:
            txt.set_fontsize(16)
        elif rid == "111":
            txt.set_fontsize(13)
        else:
            txt.set_fontsize(14)
        txt.set_fontweight("bold")
        txt.set_color(label_color[rid])
        x, y = txt.get_position()
        dx, dy = nudges.get(rid, (0.0, 0.0))
        txt.set_position((x + dx, y + dy))

    ax.set_aspect("equal", adjustable="box")
    ax.set_axis_off()

    fig.text(0.5, 0.96, args.title, ha="center", va="center", fontsize=17, fontweight="bold", color="#222222")
    if args.subtitle:
        fig.text(0.5, 0.915, args.subtitle, ha="center", va="center", fontsize=9.5,
                 fontweight="bold", color="#333333")

    # External labels placed like the Python 25e design (left / right / lower-right).
    set_lab_fs = 14.5
    fig.text(0.08, 0.85, f"{args.a_name}\n({args.a_total:,})", ha="left", va="center",
             fontsize=set_lab_fs, fontweight="bold", color="#2A2A2A")
    fig.text(0.76, 0.85, f"{args.b_name}\n({args.b_total:,})", ha="left", va="center",
             fontsize=set_lab_fs, fontweight="bold", color="#2A2A2A")
    fig.text(0.74, 0.14, f"{args.c_name}\n({args.c_total:,})", ha="left", va="center",
             fontsize=set_lab_fs, fontweight="bold", color="#2A2A2A")

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=300, bbox_inches="tight", facecolor=fig.get_facecolor())
    plt.close(fig)


if __name__ == "__main__":
    main()
