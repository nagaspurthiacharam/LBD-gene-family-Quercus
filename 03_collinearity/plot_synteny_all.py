#!/usr/bin/env python3
"""
Unified Synteny Plotter — LBD Gene Family in Quercus rubra
==========================================================
Generates publication-quality dual-track synteny plots for FOUR comparisons:

    1. Q. rubra  vs  Q. suber           (QrQs)   LG  vs  NW
    2. Q. rubra  vs  A. thaliana        (QrAt)   LG  vs  CP
    3. Q. rubra  vs  P. trichocarpa     (QrPt)   LG  vs  NC
    4. Q. rubra  vs  O. sativa          (QrOs)   LG  vs  CM

Highlights the 42 QrLBD genes (and their Q. suber XP_ orthologues) in red.
Scientific names are typeset in 12-pt italic.

Folder layout expected (script reads from its own directory):

    Synteny/
        plot_synteny_all.py            <-- this file
        LBD_highlight_QrQs.txt         <-- 42 QrLBD ids + 37 QsLBD XP_ accns
        Qr_Qs/    QrQs.gff   QrQs.collinearity   QrQs.ctl
        At_Qr/    QrAt.gff   QrAt.collinearity   QrAt.ctl
        Qr_Pt/    QrPt.gff   QrPt.collinearity   QrPt.ctl
        Qr_os/    QrOs.gff   QrOs.collinearity   QrOs.ctl
        plots/    <-- output PNG / PDF

Run:    python plot_synteny_all.py
Author: Worlasi (rewritten for the 42-gene QrLBD set, 2026)
"""

import os, re, sys
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.path as mpath


# =============================================================
#  CONFIG
# =============================================================

# BASE = folder that contains this script
BASE    = os.path.dirname(os.path.abspath(__file__))
OUT_DIR = os.path.join(BASE, "plots")
os.makedirs(OUT_DIR, exist_ok=True)

LBD_FILE = os.path.join(BASE, "LBD_highlight_QrQs.txt")

# species-pair definitions
COMPARISONS = {
    "QrQs": {
        "subdir":     "Qr_Qs",
        "stem":       "QrQs",
        "sp1_prefix": "LG",  "sp2_prefix": "NW",
        "sp1_name":   "Quercus rubra",
        "sp2_name":   "Quercus suber",
        "sp2_color":  "#E8913A",
        "sp2_edge":   "#B5691E",
        "tight":      True,   # many scaffolds → no gaps
    },
    "QrAt": {
        "subdir":     "At_Qr",
        "stem":       "QrAt",
        "sp1_prefix": "LG",  "sp2_prefix": "CP",
        "sp1_name":   "Quercus rubra",
        "sp2_name":   "Arabidopsis thaliana",
        "sp2_color":  "#7CB342",
        "sp2_edge":   "#33691E",
        "tight":      False,
    },
    "QrPt": {
        "subdir":     "Qr_Pt",
        "stem":       "QrPt",
        "sp1_prefix": "LG",  "sp2_prefix": "NC",
        "sp1_name":   "Quercus rubra",
        "sp2_name":   "Populus trichocarpa",
        "sp2_color":  "#9575CD",
        "sp2_edge":   "#4527A0",
        "tight":      False,
    },
    "QrOs": {
        "subdir":     "Qr_os",
        "stem":       "QrOs",
        "sp1_prefix": "LG",  "sp2_prefix": "CM",
        "sp1_name":   "Quercus rubra",
        "sp2_name":   "Oryza sativa",
        "sp2_color":  "#FFB74D",
        "sp2_edge":   "#E65100",
        "tight":      False,
    },
}


# =============================================================
#  HELPERS
# =============================================================

def natural_sort_key(s):
    return [int(x) if x.isdigit() else x.lower()
            for x in re.split(r"(\d+)", s)]


def parse_gff(path):
    """Returns gene_info {id: (chrom,start,end)} and chrom_max {chrom: max_end}."""
    gene_info, chrom_max = {}, defaultdict(int)
    with open(path) as fh:
        for line in fh:
            p = line.strip().split("\t")
            if len(p) < 4:
                continue
            chrom, gid = p[0], p[1]
            try:
                s, e = int(p[2]), int(p[3])
            except ValueError:
                continue
            gene_info[gid] = (chrom, s, e)
            if e > chrom_max[chrom]:
                chrom_max[chrom] = e
    return gene_info, chrom_max


def parse_collinearity(path, p1, p2):
    """Return list of inter-species blocks; chr1 is always sp1, chr2 is sp2."""
    blocks, current = [], None
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith("##"):
                m = re.search(r"(\S+)&(\S+)\s+(plus|minus)", line)
                current = None
                if not m:
                    continue
                c1, c2, ori = m.group(1), m.group(2), m.group(3)
                if c1.startswith(p1) and c2.startswith(p2):
                    pass
                elif c1.startswith(p2) and c2.startswith(p1):
                    c1, c2 = c2, c1
                else:
                    continue
                current = {"chr1": c1, "chr2": c2,
                           "pairs": [], "orientation": ori}
                blocks.append(current)
            elif current is not None and not line.startswith("#"):
                p = line.split("\t")
                if len(p) >= 3:
                    current["pairs"].append((p[1].strip(), p[2].strip()))
    return blocks


def load_lbd_genes(path):
    if not os.path.exists(path):
        return set()
    with open(path) as fh:
        return {ln.strip() for ln in fh if ln.strip()}


def shorten(chrom, prefix):
    """Compact chromosome / scaffold display label."""
    if prefix == "LG":
        return chrom                              # LG1 ... LG12
    if prefix == "CP":                            # Arabidopsis
        cp = {"CP002684_1": "Chr1", "CP002685_1": "Chr2",
              "CP002686_1": "Chr3", "CP002687_1": "Chr4",
              "CP002688_1": "Chr5"}
        return cp.get(chrom, chrom)
    if prefix == "NC":                            # Populus
        m = re.search(r"NC_(\d+)", chrom)
        if m:
            n = int(m.group(1))
            if n == 9143:
                return "Mt"
            return f"Chr{n - 37284:02d}"          # NC_037285 → Chr1
        return chrom
    if prefix == "NW":                            # Q. suber scaffolds
        m = re.search(r"NW_(\d+)", chrom)
        return f"S{int(m.group(1)) - 27069000}" if m else chrom
    if prefix == "CM":                            # Rice
        m = re.search(r"CM(\d+)", chrom)
        if m:
            n = int(m.group(1))
            return f"Chr{n - 125}"                # CM000126 → Chr1
        return chrom
    return chrom


# =============================================================
#  PLOT ONE COMPARISON
# =============================================================

def plot_one(name, cfg, lbd_genes):
    print(f"\n{'='*60}")
    print(f"  {name}: {cfg['sp1_name']}  vs  {cfg['sp2_name']}")
    print(f"{'='*60}")

    gff_path  = os.path.join(BASE, cfg["subdir"], cfg["stem"] + ".gff")
    coll_path = os.path.join(BASE, cfg["subdir"], cfg["stem"] + ".collinearity")
    if not (os.path.exists(gff_path) and os.path.exists(coll_path)):
        print(f"  SKIP: missing input file(s) under {cfg['subdir']}/")
        return

    p1, p2 = cfg["sp1_prefix"], cfg["sp2_prefix"]

    gene_info, chrom_max = parse_gff(gff_path)
    print(f"  GFF: {len(gene_info)} genes, {len(chrom_max)} sequences")

    blocks = parse_collinearity(coll_path, p1, p2)
    blocks = [b for b in blocks if b["pairs"]]
    total_pairs = sum(len(b["pairs"]) for b in blocks)
    print(f"  Inter-species blocks: {len(blocks)}  |  pairs: {total_pairs}")

    if not blocks:
        print("  No inter-species blocks; skipping.")
        return

    # only keep sp2 sequences that actually carry blocks
    sp1_chroms = sorted([c for c in chrom_max if c.startswith(p1)],
                        key=natural_sort_key)
    sp2_with_blocks = sorted({b["chr2"] for b in blocks},
                             key=natural_sort_key)
    print(f"  sp1: {len(sp1_chroms)}  |  sp2 (with blocks): {len(sp2_with_blocks)}")

    sp1_sizes = [chrom_max[c] for c in sp1_chroms]
    sp2_sizes = [chrom_max[c] for c in sp2_with_blocks]
    total_sp1 = sum(sp1_sizes) or 1
    total_sp2 = sum(sp2_sizes) or 1

    # ------- figure dimensions -------
    if cfg["tight"] or len(sp2_with_blocks) > 100:
        plot_width = 30; gap_frac = 0.0
    elif len(sp2_with_blocks) > 30:
        plot_width = 24; gap_frac = 0.001
    else:
        plot_width = 18; gap_frac = 0.012

    chrom_h, link_space = 0.5, 7.5
    margin_top, margin_bot = 2.0, 2.0
    plot_height = margin_top + chrom_h + link_space + chrom_h + margin_bot

    fig, ax = plt.subplots(1, 1, figsize=(plot_width, plot_height))

    y_sp1 = margin_top
    y_sp2 = margin_top + chrom_h + link_space

    # ------- chromosome layout -------
    def layout(chroms, sizes, total):
        n = len(chroms)
        gap = gap_frac * plot_width
        avail = plot_width - gap * max(n - 1, 0)
        pos = {}; x = 0.0
        for i, (c, sz) in enumerate(zip(chroms, sizes)):
            min_w = 0.01 if n > 100 else 0.05
            w = max((sz / total) * avail, min_w)
            pos[c] = (x, x + w)
            x += w + (gap if i < n - 1 else 0)
        return pos

    sp1_pos = layout(sp1_chroms,      sp1_sizes, total_sp1)
    sp2_pos = layout(sp2_with_blocks, sp2_sizes, total_sp2)

    # ------- chromosome rectangles -------
    def draw_track(chroms, pos, y, on_top, face, edge):
        many = len(chroms)
        for c in chroms:
            xs, xe = pos[c]
            w = xe - xs
            lw = 0.3 if many > 100 else (0.6 if many > 30 else 1.0)
            rect = patches.FancyBboxPatch(
                (xs, y), w, chrom_h,
                boxstyle="round,pad=0.01",
                facecolor=face, edgecolor=edge, linewidth=lw, alpha=0.85)
            ax.add_patch(rect)

            short = shorten(c, p1 if on_top else p2)
            if many > 100:   fs, rot = 2.0, 90
            elif many > 30:  fs, rot = 4.0, 90
            elif many > 15:  fs, rot = 6.0, 45
            else:            fs, rot = 8.5, 0

            if on_top:
                ax.text(xs + w/2, y - 0.15, short,
                        ha="center", va="top", fontsize=fs,
                        rotation=rot, fontweight="bold")
            else:
                ax.text(xs + w/2, y + chrom_h + 0.15, short,
                        ha="center", va="bottom", fontsize=fs,
                        rotation=rot, fontweight="bold")

    draw_track(sp1_chroms,      sp1_pos, y_sp1, True,  "#4A90D9", "#2C5F8A")
    draw_track(sp2_with_blocks, sp2_pos, y_sp2, False, cfg["sp2_color"], cfg["sp2_edge"])

    # ------- gene → x mapping -------
    def gene_x(gid):
        if gid not in gene_info:
            return None, None
        ch, s, e = gene_info[gid]
        mid = 0.5 * (s + e)
        mx  = chrom_max.get(ch, 1) or 1
        frac = mid / mx
        if ch in sp1_pos:
            xs, xe = sp1_pos[ch]; return xs + frac * (xe - xs), "sp1"
        if ch in sp2_pos:
            xs, xe = sp2_pos[ch]; return xs + frac * (xe - xs), "sp2"
        return None, None

    # ------- adaptive transparency -------
    if   total_pairs > 5000: base_a, base_lw = 0.05, 0.25
    elif total_pairs > 2000: base_a, base_lw = 0.08, 0.40
    elif total_pairs >  500: base_a, base_lw = 0.15, 0.50
    else:                    base_a, base_lw = 0.30, 0.60

    # ------- draw bezier ribbons -------
    y1_link = y_sp1 + chrom_h         # bottom of top track
    y2_link = y_sp2                   # top of bottom track
    mid_y   = 0.5 * (y1_link + y2_link)

    lbd_links = []; drawn = 0; skipped = 0

    for blk in blocks:
        for g1, g2 in blk["pairs"]:
            x1, s1 = gene_x(g1)
            x2, s2 = gene_x(g2)
            if x1 is None or x2 is None:
                skipped += 1; continue
            if   s1 == "sp2" and s2 == "sp1": x1, x2 = x2, x1
            elif s1 == s2:
                skipped += 1; continue

            if g1 in lbd_genes or g2 in lbd_genes:
                lbd_links.append((x1, x2))
                continue

            verts = [(x1, y1_link),
                     (x1, mid_y - 0.4),
                     (x2, mid_y + 0.4),
                     (x2, y2_link)]
            path = mpath.Path(verts, [1, 4, 4, 4])
            ax.add_patch(patches.PathPatch(
                path, facecolor="none", edgecolor="#888888",
                alpha=base_a, linewidth=base_lw))
            drawn += 1

    # LBD ribbons drawn last → on top
    for x1, x2 in lbd_links:
        verts = [(x1, y1_link),
                 (x1, mid_y - 0.4),
                 (x2, mid_y + 0.4),
                 (x2, y2_link)]
        path = mpath.Path(verts, [1, 4, 4, 4])
        ax.add_patch(patches.PathPatch(
            path, facecolor="none", edgecolor="#E74C3C",
            alpha=0.95, linewidth=1.8))
        drawn += 1

    print(f"  ribbons drawn: {drawn}  (skipped {skipped})")
    print(f"  LBD ribbons (red): {len(lbd_links)}")

    # ------- title (LARGE italic species names) -------
    sp1_tex = r"$\mathit{" + cfg["sp1_name"].replace(" ", r"\ ") + "}$"
    sp2_tex = r"$\mathit{" + cfg["sp2_name"].replace(" ", r"\ ") + "}$"
    title_str = f"{sp1_tex}   vs   {sp2_tex}"

    ax.text(plot_width / 2, 0.45, title_str,
            ha="center", va="center", fontsize=28, fontweight="bold")

    subtitle = (f"{len(blocks)} collinear blocks  |  {total_pairs} gene pairs"
                f"  |  {len(sp2_with_blocks)} {p2} sequence(s)")
    if lbd_links:
        subtitle += f"  |  {len(lbd_links)} LBD pair(s) in red"
    ax.text(plot_width / 2, 1.30, subtitle,
            ha="center", va="center", fontsize=11, color="#555555")

    # ------- side labels (LARGE italic, prominent) -------
    ax.text(-0.6, y_sp1 + chrom_h / 2, cfg["sp1_name"],
            ha="right", va="center",
            fontsize=22, style="italic", fontweight="bold",
            color="#2C5F8A")
    ax.text(-0.6, y_sp2 + chrom_h / 2, cfg["sp2_name"],
            ha="right", va="center",
            fontsize=22, style="italic", fontweight="bold",
            color=cfg["sp2_edge"])

    # ------- finalise -------
    # extra left margin so the large italic species name labels don't clip
    ax.set_xlim(-6.5, plot_width + 0.5)
    ax.set_ylim(plot_height + 0.4, -0.2)
    ax.set_aspect("auto")
    ax.axis("off")

    out_png = os.path.join(OUT_DIR, f"{name}_synteny.png")
    out_pdf = os.path.join(OUT_DIR, f"{name}_synteny.pdf")
    fig.savefig(out_png, dpi=300, bbox_inches="tight", facecolor="white")
    fig.savefig(out_pdf,            bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"  saved: {out_png}")
    print(f"  saved: {out_pdf}")


# =============================================================
#  MAIN
# =============================================================

def main():
    print("=" * 60)
    print("  Synteny Plotter — LBD Gene Family in Quercus rubra")
    print("  42 QrLBD genes (incl. 1 sub-isoform) + 37 QsLBD orthologues")
    print("=" * 60)
    print(f"  Working in:  {BASE}")
    print(f"  Output dir:  {OUT_DIR}")

    lbd = load_lbd_genes(LBD_FILE)
    n_qr = sum(1 for g in lbd if g.startswith("Qurub"))
    n_qs = sum(1 for g in lbd if g.startswith("XP_"))
    print(f"  LBD highlight list: {len(lbd)} ids "
          f"({n_qr} QrLBD + {n_qs} QsLBD orthologues)")

    for key, cfg in COMPARISONS.items():
        plot_one(key, cfg, lbd)

    print("\n" + "=" * 60)
    print(f"  All four plots written to {OUT_DIR}")
    print("=" * 60)


if __name__ == "__main__":
    main()
