#!/usr/bin/env python3
"""Annotate peak midpoints and render genomic-distribution pies."""

from __future__ import annotations

from pathlib import Path
import re
import subprocess
import sys
from typing import Dict, List, Tuple, Optional, Set

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

CUTRUN_ROOT_CANDIDATES = [Path.cwd() / "cutrun_work"]

ASSAYS = ["MCM3", "NONO", "PSPC1"]

USE_NT_KD_UNION = True
PEAK_UNIT_LABEL = "peaks"

BEDTOOLS = Path("/opt/anaconda3/envs/cutrun_env/bin/bedtools")

# Peak distribution annotations
PROMOTER_UP_BP = 1000
PROMOTER_DOWN_BP = 1000
DOWNSTREAM_BP = 300

GTF_CANDIDATES = [
    Path("ref/gencode.vM25.annotation.gtf"),
    Path("ref/gencode.vM25.annotation.gtf.gz"),
    Path("ref/gencode.vM23.annotation.gtf"),
    Path("ref/gencode.vM23.annotation.gtf.gz"),
    Path("ref/gencode.annotation.gtf"),
    Path("ref/gencode.annotation.gtf.gz"),
]

PIE_EDGE_COLOR = "black"
PIE_EDGE_LW = 1.6

COLORS = {
    "Promoter (±1kb)": "#1f78b4",
    "5' UTR": "#33a02c",
    "3' UTR": "#fb9a99",
    "1st Exon": "#e31a1c",
    "Other Exon": "#fdbf6f",
    "1st Intron": "#cab2d6",
    "Other Intron": "#6a3d9a",
    "Downstream (<=300bp)": "#ffff99",
    "Distal Intergenic": "#b15928",
}

CATEGORY_ORDER = [
    "Promoter (±1kb)",
    "5' UTR",
    "3' UTR",
    "1st Exon",
    "Other Exon",
    "1st Intron",
    "Other Intron",
    "Downstream (<=300bp)",
    "Distal Intergenic",
]

def die(msg: str, code: int = 1) -> None:
    print(f"[ERROR] {msg}", file=sys.stderr)
    sys.exit(code)

def run(cmd: str) -> None:
    p = subprocess.run(cmd, shell=True, executable="/bin/bash")
    if p.returncode != 0:
        die(f"command failed: {cmd}")

def first_existing_dir(paths: List[Path]) -> Optional[Path]:
    for p in paths:
        if p.exists() and p.is_dir():
            return p.resolve()
    return None

WORK_ROOT = CUTRUN_ROOT_CANDIDATES[0]

OUTDIR = WORK_ROOT / "data" / "figure_inputs" / "peak_distribution" / "intermediate"
TMPDIR = OUTDIR / "intermediate"
PAPER_TMPDIR = OUTDIR / "paperfig_intermediate"

CONSENSUS_BASE = WORK_ROOT / "05_consensus_macs2"
BEDGRAPH_ROOT = WORK_ROOT / "03_bedgraph"


def ensure_dirs() -> None:
    if not BEDTOOLS.exists():
        die(f"bedtools not found at: {BEDTOOLS}")
    OUTDIR.mkdir(parents=True, exist_ok=True)
    TMPDIR.mkdir(parents=True, exist_ok=True)
    PAPER_TMPDIR.mkdir(parents=True, exist_ok=True)

def sort_bed_inplace(bed: Path) -> None:
    if (not bed.exists()) or bed.stat().st_size == 0:
        return
    tmp = bed.with_suffix(bed.suffix + ".sorted.tmp")
    run(f'sort -k1,1 -k2,2n "{bed}" > "{tmp}" && mv "{tmp}" "{bed}"')

def autodetect_gtf() -> Path:
    for p in GTF_CANDIDATES:
        if p.exists():
            print(f"[OK] Using GTF: {p}")
            return p
    die("Could not autodetect GTF. Tried:\n  - " + "\n  - ".join(str(p) for p in GTF_CANDIDATES))

def _sanitize_intervals_df(df: pd.DataFrame, chrom="chrom", start="start", end="end") -> pd.DataFrame:
    if df.empty:
        return df.copy()
    out = df.copy()
    out[chrom] = out[chrom].astype(str)
    out = out[out[chrom].str.len() > 0].copy()

    out[start] = pd.to_numeric(out[start], errors="coerce").fillna(0).astype(np.int64)
    out[end] = pd.to_numeric(out[end], errors="coerce").fillna(0).astype(np.int64)

    out.loc[out[start] < 0, start] = 0
    out.loc[out[end] < 0, end] = 0

    bad = out[end] <= out[start]
    out.loc[bad, end] = out.loc[bad, start] + 1
    out = out[out[end] > out[start]].copy()
    return out

def write_bed3(df: pd.DataFrame, out_bed: Path) -> None:
    if df.empty:
        out_bed.write_text("")
        return
    d = _sanitize_intervals_df(df[["chrom", "start", "end"]].copy())
    if d.empty:
        out_bed.write_text("")
        return
    d.to_csv(out_bed, sep="\t", header=False, index=False)
    sort_bed_inplace(out_bed)

def count_bed_lines(p: Path) -> int:
    if (not p.exists()) or p.stat().st_size == 0:
        return 0
    with open(p, "r") as f:
        return sum(1 for line in f if line.strip() and not line.startswith("#"))

def read_bed3_as_set(p: Path) -> Set[Tuple[str, int, int]]:
    if (not p.exists()) or p.stat().st_size == 0:
        return set()
    out: Set[Tuple[str, int, int]] = set()
    with open(p, "r") as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            chrom = str(parts[0])
            try:
                start = int(float(parts[1]))
                end = int(float(parts[2]))
            except Exception:
                continue
            if end <= start:
                continue
            out.add((chrom, start, end))
    return out


# =============================================================================
# PEAK FINDERS (authoritative peak-set definition)
# =============================================================================
def _find_nt_bed(assay: str) -> Path:
    p = CONSENSUS_BASE / "NT" / assay / "union.merge.bed"
    if p.exists():
        return p
    # fallback search
    beds = sorted(CONSENSUS_BASE.rglob("union.merge.bed"))
    beds = [b for b in beds if f"/NT/{assay}/union.merge.bed" in str(b).replace("\\", "/")]
    if not beds:
        die(f"No NT union.merge.bed found for assay={assay} under: {CONSENSUS_BASE}")
    beds.sort(key=lambda x: (len(str(x)), str(x)))
    return beds[0]

def _find_kd_bed(assay: str) -> Optional[Path]:
    cands = [
        CONSENSUS_BASE / f"sh{assay}" / assay / "union.merge.bed",
        CONSENSUS_BASE / "KD" / assay / "union.merge.bed",
        CONSENSUS_BASE / f"sh{assay}" / assay / "consensus.min1.bed",
    ]
    for c in cands:
        if c.exists():
            return c
    return None

def build_union_bed_unique_rows(nt: Path, kd: Path, out_union: Path) -> Path:
    df_nt = pd.read_csv(nt, sep="\t", header=None, usecols=[0, 1, 2], names=["chrom", "start", "end"])
    df_kd = pd.read_csv(kd, sep="\t", header=None, usecols=[0, 1, 2], names=["chrom", "start", "end"])

    df = pd.concat([df_nt, df_kd], ignore_index=True)
    df = _sanitize_intervals_df(df)
    if df.empty:
        out_union.write_text("")
        return out_union

    df = df.drop_duplicates(subset=["chrom", "start", "end"]).copy()
    df.to_csv(out_union, sep="\t", header=False, index=False)
    sort_bed_inplace(out_union)
    return out_union

def find_peaks_bed(assay: str) -> Path:
    nt = _find_nt_bed(assay)
    if not USE_NT_KD_UNION:
        return nt
    kd = _find_kd_bed(assay)
    if kd is None:
        print(f"[WARN] KD bed not found for {assay}; using NT only.")
        return nt
    out_union = TMPDIR / f"{assay}.NT_union_KD.unique_rows.bed"
    build_union_bed_unique_rows(nt, kd, out_union)
    if (not out_union.exists()) or out_union.stat().st_size == 0:
        die(f"Union bed produced empty file for {assay}: {out_union}")
    print(f"[OK] Using NT∪KD unique-row union for {assay}:\n     NT={nt}\n     KD={kd}\n     UNION={out_union}")
    return out_union


# =============================================================================
# GTF parsing & canonical transcript (PeakDist + Paperfigs shared)
# =============================================================================
def gtf_parse(gtf: Path) -> pd.DataFrame:
    cols = ["chrom", "source", "feature", "start", "end", "score", "strand", "frame", "attr"]
    df = pd.read_csv(
        gtf, sep="\t", comment="#", header=None, names=cols,
        dtype={"chrom": str, "feature": str, "strand": str, "attr": str},
        low_memory=False,
    )
    df["start"] = pd.to_numeric(df["start"], errors="coerce").fillna(1).astype(np.int64) - 1
    df["end"] = pd.to_numeric(df["end"], errors="coerce").fillna(1).astype(np.int64)

    def get_attr(s: str, key: str) -> str:
        m = re.search(rf'{key} "([^"]+)"', str(s))
        return m.group(1) if m else ""

    df["gene_id"] = df["attr"].map(lambda x: get_attr(x, "gene_id"))
    df["transcript_id"] = df["attr"].map(lambda x: get_attr(x, "transcript_id"))
    return df

def harmonize_chrom_style_from_bed(peaks_bed: Path, df_gtf: pd.DataFrame) -> pd.DataFrame:
    # Infer naming from canonical chromosomes, not from the first record. The
    # controlled peak sets mix unprefixed alt contigs such as GL456210.1 with
    # chr-prefixed canonical chromosomes.
    has_chr: Optional[bool] = None
    with open(peaks_bed, "r") as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            chrom = line.split("\t", 1)[0]
            if re.fullmatch(r"chr(?:[1-9]|1[0-9]|X|Y|M)", chrom):
                has_chr = True
                break
            if re.fullmatch(r"(?:[1-9]|1[0-9]|X|Y|M|MT)", chrom):
                has_chr = False
                break
    if has_chr is None:
        die(f"Could not infer canonical chromosome style from: {peaks_bed}")
    out = df_gtf.copy()
    if has_chr:
        out["chrom"] = out["chrom"].astype(str).apply(lambda x: x if x.startswith("chr") else "chr" + x)
    else:
        out["chrom"] = out["chrom"].astype(str).str.replace(r"^chr", "", regex=True)
    return out

def pick_canonical_transcripts(df: pd.DataFrame) -> pd.DataFrame:
    ex = df[df["feature"] == "exon"].copy()
    ex = ex[(ex["gene_id"] != "") & (ex["transcript_id"] != "")]
    if ex.empty:
        die("No exon features with transcript_id found in GTF.")
    ex["len"] = (ex["end"] - ex["start"]).clip(lower=0)
    sums = ex.groupby(["gene_id", "transcript_id"], as_index=False)["len"].sum()
    sums = sums.sort_values(["gene_id", "len", "transcript_id"], ascending=[True, False, True])
    return sums.drop_duplicates("gene_id", keep="first")[["gene_id", "transcript_id"]].copy()

def build_transcript_tss_tes(df_gtf: pd.DataFrame, canon: pd.DataFrame) -> pd.DataFrame:
    tr = df_gtf[df_gtf["feature"] == "transcript"].copy()
    tr = tr.merge(canon, on=["gene_id", "transcript_id"], how="inner")
    if tr.empty:
        ex = build_transcript_exons(df_gtf, canon)
        if ex.empty:
            die("Cannot derive transcript extents; no canonical exons found.")
        tr = ex.groupby(["chrom", "strand", "gene_id", "transcript_id"], as_index=False).agg(
            start=("start", "min"),
            end=("end", "max"),
        )
    else:
        tr = tr[["chrom", "start", "end", "strand", "gene_id", "transcript_id"]].copy()

    tr["tss"] = np.where(tr["strand"] == "+", tr["start"], tr["end"]).astype(np.int64)
    tr["tes"] = np.where(tr["strand"] == "+", tr["end"], tr["start"]).astype(np.int64)
    return tr

def build_transcript_exons(df_gtf: pd.DataFrame, canon: pd.DataFrame) -> pd.DataFrame:
    ex = df_gtf[df_gtf["feature"] == "exon"].copy()
    ex = ex.merge(canon, on=["gene_id", "transcript_id"], how="inner")
    ex = ex[["chrom", "start", "end", "strand", "gene_id", "transcript_id"]].copy()
    ex = _sanitize_intervals_df(ex)
    return ex

def build_transcript_utr(df_gtf: pd.DataFrame, canon: pd.DataFrame, tr: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
    u5 = df_gtf[df_gtf["feature"].isin(["five_prime_utr", "5UTR", "five_prime_UTR"])].copy()
    u3 = df_gtf[df_gtf["feature"].isin(["three_prime_utr", "3UTR", "three_prime_UTR"])].copy()

    if (not u5.empty) or (not u3.empty):
        u5 = _sanitize_intervals_df(u5.merge(canon, on=["gene_id", "transcript_id"], how="inner")[["chrom","start","end"]])
        u3 = _sanitize_intervals_df(u3.merge(canon, on=["gene_id", "transcript_id"], how="inner")[["chrom","start","end"]])
        return u5, u3

    utr = df_gtf[df_gtf["feature"] == "UTR"].copy()
    if utr.empty:
        return (pd.DataFrame(columns=["chrom","start","end"]), pd.DataFrame(columns=["chrom","start","end"]))

    utr = utr.merge(canon, on=["gene_id", "transcript_id"], how="inner")
    if utr.empty:
        return (pd.DataFrame(columns=["chrom","start","end"]), pd.DataFrame(columns=["chrom","start","end"]))

    utr = utr.merge(tr[["gene_id","transcript_id","tss","tes"]], on=["gene_id","transcript_id"], how="inner")
    utr = _sanitize_intervals_df(utr)
    if utr.empty:
        return (pd.DataFrame(columns=["chrom","start","end"]), pd.DataFrame(columns=["chrom","start","end"]))

    mid = ((utr["start"].astype(np.int64) + utr["end"].astype(np.int64)) // 2).astype(np.int64)
    tss = utr["tss"].astype(np.int64)
    tes = utr["tes"].astype(np.int64)

    is_u5 = (mid - tss).abs() <= (mid - tes).abs()
    u5_df = _sanitize_intervals_df(utr.loc[is_u5, ["chrom","start","end"]].copy())
    u3_df = _sanitize_intervals_df(utr.loc[~is_u5, ["chrom","start","end"]].copy())
    return u5_df, u3_df

def build_promoter_bed(tr: pd.DataFrame) -> pd.DataFrame:
    s = (tr["tss"] - PROMOTER_UP_BP).clip(lower=0).astype(np.int64)
    e = (tr["tss"] + PROMOTER_DOWN_BP).astype(np.int64)
    df = pd.DataFrame({"chrom": tr["chrom"].astype(str), "start": s, "end": e})
    return _sanitize_intervals_df(df)

def build_downstream_bed(tr: pd.DataFrame) -> pd.DataFrame:
    strand = tr["strand"].astype(str)
    tes = tr["tes"].astype(np.int64)
    s = np.where(strand == "+", tes, (tes - DOWNSTREAM_BP).clip(lower=0)).astype(np.int64)
    e = np.where(strand == "+", tes + DOWNSTREAM_BP, tes).astype(np.int64)
    df = pd.DataFrame({"chrom": tr["chrom"].astype(str), "start": s, "end": e})
    return _sanitize_intervals_df(df)

def build_first_other_exons(ex: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
    if ex.empty:
        return (pd.DataFrame(columns=["chrom","start","end"]), pd.DataFrame(columns=["chrom","start","end"]))
    ex2 = ex.copy()
    ex2["_first_key"] = np.where(ex2["strand"] == "+", ex2["start"], -ex2["end"]).astype(np.float64)
    idx = ex2.groupby(["gene_id","transcript_id"])["_first_key"].idxmin()
    first = _sanitize_intervals_df(ex2.loc[idx, ["chrom","start","end"]].copy())
    other = _sanitize_intervals_df(ex2.drop(index=idx)[["chrom","start","end"]].copy())
    return first, other

def build_introns_from_exons(ex: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
    if ex.empty:
        return (pd.DataFrame(columns=["chrom","start","end"]), pd.DataFrame(columns=["chrom","start","end"]))

    first_introns: List[Tuple[str,int,int]] = []
    other_introns: List[Tuple[str,int,int]] = []

    for (_, _), g in ex.groupby(["gene_id","transcript_id"]):
        strand = str(g["strand"].iloc[0])
        gg = g[["chrom","start","end"]].copy().sort_values(["start","end"])
        exons = list(gg.itertuples(index=False, name=None))
        if len(exons) < 2:
            continue

        intrs: List[Tuple[str,int,int]] = []
        for j in range(len(exons)-1):
            c1,s1,e1 = exons[j]
            c2,s2,e2 = exons[j+1]
            if c2 != c1:
                continue
            if int(s2) > int(e1):
                intrs.append((str(c1), int(e1), int(s2)))

        intrs = [(c,s,e) for (c,s,e) in intrs if e > s and s >= 0]
        if not intrs:
            continue

        if strand == "+":
            first_introns.append(intrs[0])
            other_introns.extend(intrs[1:])
        else:
            first_introns.append(intrs[-1])
            other_introns.extend(intrs[:-1])

    fi = _sanitize_intervals_df(pd.DataFrame(first_introns, columns=["chrom","start","end"]))
    oi = _sanitize_intervals_df(pd.DataFrame(other_introns, columns=["chrom","start","end"]))
    return fi, oi


# =============================================================================
# PEAKDIST: assignment via bedtools
# =============================================================================
def bed_to_peak_midpoints(peaks_bed: Path, out_bed: Path) -> int:
    df = pd.read_csv(peaks_bed, sep="\t", header=None, usecols=[0, 1, 2], names=["chrom", "start", "end"])
    if df.empty:
        out_bed.write_text("")
        return 0

    df = _sanitize_intervals_df(df)
    if df.empty:
        out_bed.write_text("")
        return 0

    mid = ((df["start"].astype(np.int64) + df["end"].astype(np.int64)) // 2).astype(np.int64)
    mid[mid < 0] = 0

    out = pd.DataFrame(
        {"chrom": df["chrom"].astype(str).values,
         "start": mid.values,
         "end": (mid + 1).values,
         "peak_id": np.arange(len(mid), dtype=np.int64)}
    )
    out = _sanitize_intervals_df(out)
    if out.empty:
        out_bed.write_text("")
        return 0

    out.to_csv(out_bed, sep="\t", header=False, index=False)
    sort_bed_inplace(out_bed)
    return int(len(out))

def bedtools_intersect_ids(a_bed4: Path, b_bed3: Path, out_tsv: Path) -> pd.Series:
    if (not a_bed4.exists()) or a_bed4.stat().st_size == 0:
        return pd.Series([], dtype=np.int64)
    if (not b_bed3.exists()) or b_bed3.stat().st_size == 0:
        return pd.Series([], dtype=np.int64)

    sort_bed_inplace(a_bed4)
    sort_bed_inplace(b_bed3)

    run(f'{BEDTOOLS} intersect -u -a "{a_bed4}" -b "{b_bed3}" > "{out_tsv}"')
    if (not out_tsv.exists()) or out_tsv.stat().st_size == 0:
        return pd.Series([], dtype=np.int64)

    df = pd.read_csv(out_tsv, sep="\t", header=None)
    ids = pd.to_numeric(df.iloc[:, 3], errors="coerce").dropna().astype(np.int64)
    return ids

def assign_categories(peaks_mid_bed: Path, category_beds: Dict[str, Path]) -> Dict[str, int]:
    assigned: Set[int] = set()
    counts: Dict[str, int] = {k: 0 for k in CATEGORY_ORDER}

    with open(peaks_mid_bed, "r") as f:
        total = sum(1 for line in f if line.strip())

    for cat in CATEGORY_ORDER[:-1]:
        ids = bedtools_intersect_ids(
            a_bed4=peaks_mid_bed,
            b_bed3=category_beds[cat],
            out_tsv=TMPDIR / f"_hit_{peaks_mid_bed.stem}__{re.sub(r'[^A-Za-z0-9]+','_',cat)}.tsv",
        )
        new = []
        for x in ids.to_list():
            ix = int(x)
            if ix not in assigned:
                new.append(ix)
        counts[cat] = len(new)
        assigned.update(new)

    counts["Distal Intergenic"] = max(0, total - sum(counts[c] for c in CATEGORY_ORDER[:-1]))
    return counts

def pie_one(ax, counts: Dict[str, int], title: str) -> None:
    labels = []
    sizes = []
    colors = []

    total = sum(int(counts.get(k, 0)) for k in CATEGORY_ORDER)
    for k in CATEGORY_ORDER:
        v = int(counts.get(k, 0))
        labels.append(k)
        sizes.append(v)
        colors.append(COLORS.get(k, "#cccccc"))

    explode = [0.0] * len(sizes)
    if total > 0 and max(sizes) > 0:
        explode[int(np.argmax(sizes))] = 0.06

    wedges, _ = ax.pie(
        sizes,
        startangle=0,
        colors=colors,
        explode=explode,
        wedgeprops=dict(edgecolor=PIE_EDGE_COLOR, linewidth=PIE_EDGE_LW),
    )
    ax.set_aspect("equal")
    ax.set_title(title, fontsize=16, fontweight="bold", pad=10)

    legend_lines = []
    for k, v in zip(labels, sizes):
        pct = (100.0 * v / total) if total else 0.0
        legend_lines.append(f"{k} ({pct:.2f}%)")

    legend = ax.legend(
        wedges,
        legend_lines,
        loc="center left",
        bbox_to_anchor=(1.02, 0.5),
        frameon=False,
        fontsize=11,
    )
    for text in legend.get_texts():
        text.set_fontweight("bold")

def save_three_panel(out_png: Path, results: Dict[str, Dict[str, int]], n_peaks: Dict[str, int]) -> None:
    fig = plt.figure(figsize=(10.5, 13.5))
    gs = fig.add_gridspec(3, 1, hspace=0.32)

    for i, assay in enumerate(ASSAYS):
        ax = fig.add_subplot(gs[i, 0])
        title = f"{assay} CUT&RUN {PEAK_UNIT_LABEL}\n(n={n_peaks[assay]:,})"
        pie_one(ax, results[assay], title)

    plt.tight_layout()
    plt.savefig(out_png, dpi=300)
    plt.close()

def save_single(out_png: Path, assay: str, counts: Dict[str, int], n: int) -> None:
    plt.figure(figsize=(10.5, 4.6))
    ax = plt.gca()
    pie_one(ax, counts, f"{assay} CUT&RUN {PEAK_UNIT_LABEL}\n(n={n:,})")
    plt.tight_layout()
    plt.savefig(out_png, dpi=300)
    plt.close()
