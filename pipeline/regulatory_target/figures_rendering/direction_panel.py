#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse


def fmt_int(x: int) -> str:
    return f"{int(x):,}"


def _args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    for k in [
        "out", "panel_letter", "factor", "center_col", "col_up", "col_down", "col_overlap", "stroke_col"
    ]:
        p.add_argument(f"--{k.replace('_','-')}", dest=k, required=True)
    for k in [
        "n_center", "n_up", "n_down", "n_only_center", "n_only_up", "n_only_down", "n_overlap_up", "n_overlap_down",
        "width_px", "height_px", "dpi",
        "fs_panel_letter", "fs_title", "fs_subtitle", "fs_num_center", "fs_num_overlap", "fs_bottom_lab", "fs_bottom_n"
    ]:
        p.add_argument(f"--{k.replace('_','-')}", dest=k, required=True, type=float)
    for k in ["stroke_lwd", "r_side", "r_center", "cx", "cy", "lx", "ly", "rx", "ry"]:
        p.add_argument(f"--{k.replace('_','-')}", dest=k, required=True, type=float)
    p.add_argument("--hide-numbers", action="store_true")
    return p.parse_args()


def _circle_patch_axes(xy, r, xy_scale, **kwargs):
    # Compensate for non-square panel aspect so displayed shapes remain circular.
    return Ellipse(xy, width=2 * r * xy_scale, height=2 * r, transform=kwargs.pop("transform"), **kwargs)


def add_overlap_lens(ax, c1xy, r1, c2xy, r2, fill, xy_scale, alpha=0.95, z=5):
    # Draw exact intersection by clipping a circle patch with the other circle path.
    base = _circle_patch_axes(c1xy, r1, xy_scale, facecolor=fill, edgecolor="none",
                              alpha=alpha, transform=ax.transAxes, zorder=z)
    clip = _circle_patch_axes(c2xy, r2, xy_scale, facecolor="none", edgecolor="none",
                              transform=ax.transAxes)
    base.set_clip_path(clip)
    ax.add_patch(base)


def main() -> None:
    a = _args()

    fig_w = a.width_px / a.dpi
    fig_h = a.height_px / a.dpi
    fig = plt.figure(figsize=(fig_w, fig_h), dpi=int(a.dpi), facecolor="white")
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_axis_off()
    xy_scale = fig_h / fig_w

    # Circle fills (opaque)
    up_patch = _circle_patch_axes((a.lx, a.ly), a.r_side, xy_scale, facecolor=a.col_up, edgecolor="none",
                                  transform=ax.transAxes, zorder=2)
    dn_patch = _circle_patch_axes((a.rx, a.ry), a.r_side, xy_scale, facecolor=a.col_down, edgecolor="none",
                                  transform=ax.transAxes, zorder=2)
    ct_patch = _circle_patch_axes((a.cx, a.cy), a.r_center, xy_scale, facecolor=a.center_col, edgecolor="none",
                                  transform=ax.transAxes, zorder=3)
    ax.add_patch(up_patch)
    ax.add_patch(dn_patch)
    ax.add_patch(ct_patch)

    # Exact overlap lenses only
    add_overlap_lens(ax, (a.lx, a.ly), a.r_side, (a.cx, a.cy), a.r_center, a.col_overlap, xy_scale, z=4)
    add_overlap_lens(ax, (a.cx, a.cy), a.r_center, (a.lx, a.ly), a.r_side, a.col_overlap, xy_scale, z=4)
    add_overlap_lens(ax, (a.rx, a.ry), a.r_side, (a.cx, a.cy), a.r_center, a.col_overlap, xy_scale, z=4)
    add_overlap_lens(ax, (a.cx, a.cy), a.r_center, (a.rx, a.ry), a.r_side, a.col_overlap, xy_scale, z=4)

    # Circle outlines on top
    for (x, y, r) in [(a.lx, a.ly, a.r_side), (a.rx, a.ry, a.r_side), (a.cx, a.cy, a.r_center)]:
        ax.add_patch(_circle_patch_axes((x, y), r, xy_scale, facecolor="none", edgecolor=a.stroke_col,
                                        linewidth=a.stroke_lwd, transform=ax.transAxes, zorder=10))

    # Titles
    ax.text(0.06, 0.93, a.panel_letter, transform=ax.transAxes, ha="center", va="center",
            fontsize=a.fs_panel_letter, fontweight="bold", color="black", zorder=30)
    ax.text(0.50, 0.94, f"{a.factor} CUT&RUN", transform=ax.transAxes, ha="center", va="center",
            fontsize=a.fs_title, fontweight="bold", color="black", zorder=30)
    if not a.hide_numbers:
        ax.text(0.50, 0.885, f"({fmt_int(a.n_center)})", transform=ax.transAxes, ha="center", va="center",
                fontsize=a.fs_subtitle, fontweight="bold", color="black", zorder=30)

    # Numbers
    if not a.hide_numbers:
        up_ox = (a.lx + a.cx) / 2 - 0.014
        up_oy = (a.ly + a.cy) / 2
        dn_ox = (a.rx + a.cx) / 2 + 0.014
        dn_oy = (a.ry + a.cy) / 2
        ax.text(a.lx - 0.03, a.ly, fmt_int(a.n_only_up), transform=ax.transAxes, ha="center", va="center",
                fontsize=a.fs_num_center, fontweight="bold", color="white", zorder=30)
        ax.text(a.cx, a.cy, fmt_int(a.n_only_center), transform=ax.transAxes, ha="center", va="center",
                fontsize=a.fs_num_center, fontweight="bold", color=a.stroke_col, zorder=30)
        ax.text(a.rx + 0.03, a.ry, fmt_int(a.n_only_down), transform=ax.transAxes, ha="center", va="center",
                fontsize=a.fs_num_center, fontweight="bold", color="white", zorder=30)
        if a.n_overlap_up > 0:
            ax.text(up_ox, up_oy, fmt_int(a.n_overlap_up), transform=ax.transAxes, ha="center", va="center",
                    fontsize=a.fs_num_overlap, fontweight="bold", color="black", zorder=30)
        if a.n_overlap_down > 0:
            ax.text(dn_ox, dn_oy, fmt_int(a.n_overlap_down), transform=ax.transAxes, ha="center", va="center",
                    fontsize=a.fs_num_overlap, fontweight="bold", color="black", zorder=30)

    # Bottom labels
    ax.text(a.lx, 0.155, f"{a.factor}-KD Up", transform=ax.transAxes, ha="center", va="center",
            fontsize=a.fs_bottom_lab, fontweight="bold", color="black", zorder=30)
    ax.text(a.rx, 0.155, f"{a.factor}-KD Down", transform=ax.transAxes, ha="center", va="center",
            fontsize=a.fs_bottom_lab, fontweight="bold", color="black", zorder=30)
    if not a.hide_numbers:
        ax.text(a.lx, 0.110, f"({fmt_int(a.n_up)})", transform=ax.transAxes, ha="center", va="center",
                fontsize=a.fs_bottom_n, fontweight="bold", color="black", zorder=30)
        ax.text(a.rx, 0.110, f"({fmt_int(a.n_down)})", transform=ax.transAxes, ha="center", va="center",
                fontsize=a.fs_bottom_n, fontweight="bold", color="black", zorder=30)

    out = Path(a.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=int(a.dpi), facecolor="white")
    plt.close(fig)


if __name__ == "__main__":
    main()
