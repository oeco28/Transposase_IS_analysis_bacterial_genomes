#!/usr/bin/env python3
"""
05_figures.py  (v2)

Publication-quality figures from the pipeline outputs. Each panel is
independent: a figure whose input is missing is skipped with a note.

v2 design changes:
  * no insets overlaying data (the old fig6 inset sat on the ECDF curves)
  * panel letters, consistent typography, constrained layout
  * a colourblind-safe palette with one accent for the focal group
  * every panel states its n; effect sizes shown with intervals, not just
    points
  * PDF (vector, for the journal) and PNG at 400 dpi (for drafts)

Usage:
  python 05_figures.py --focal 12 --isout is_out --meta metadata.tsv \\
      --out figures
  python 05_figures.py --focal 12 --only fig6
"""

import argparse
import os

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

# Okabe-Ito derived: safe for the common colour vision deficiencies, and
# still legible in greyscale because the accent is much darker than the
# neutral. Colour carries the focal/other contrast and nothing else.
FOCAL = "#B0353C"
FOCAL_L = "#E8B4B7"
NEUTRAL = "#7A7A7A"
NEUTRAL_L = "#D6D6D6"
INK = "#1A1A1A"
MUTED = "#6E6E6E"

plt.rcParams.update({
    "figure.dpi": 110,
    "savefig.dpi": 400,
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.03,
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
    "font.size": 8.5,
    "axes.titlesize": 9.5,
    "axes.labelsize": 8.5,
    "axes.labelcolor": INK,
    "axes.edgecolor": "#4D4D4D",
    "axes.linewidth": 0.7,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.titlelocation": "left",
    "axes.titlepad": 6,
    "xtick.color": MUTED,
    "ytick.color": MUTED,
    "xtick.labelsize": 7.5,
    "ytick.labelsize": 7.5,
    "xtick.major.width": 0.7,
    "ytick.major.width": 0.7,
    "xtick.major.size": 3,
    "ytick.major.size": 3,
    "legend.frameon": False,
    "legend.fontsize": 7.5,
    "legend.handlelength": 1.4,
    "figure.facecolor": "white",
    "axes.grid": False,
})


def direction_note(ax, left_text, right_text, y=-0.20):
    """
    Label which end of the axis means "more IS-associated".
    fig5 (distance) and fig6 (proximity fraction) run in OPPOSITE
    directions, so an unlabelled reader will conclude the two figures
    contradict each other when they agree.
    """
    ax.annotate(left_text, xy=(0, y), xycoords="axes fraction", fontsize=6.5,
                color=MUTED, ha="left", va="top", annotation_clip=False)
    ax.annotate(right_text, xy=(1, y), xycoords="axes fraction", fontsize=6.5,
                color=MUTED, ha="right", va="top", annotation_clip=False)


def panel_letter(ax, letter):
    ax.annotate(letter, xy=(0, 1), xycoords="axes fraction",
                xytext=(-26, 10), textcoords="offset points",
                fontsize=11, fontweight="bold", color=INK, va="top")


def save(fig, out, name):
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(out, f"{name}.{ext}"))
    plt.close(fig)
    print(f"  wrote {name}.pdf / .png")


def need(*paths):
    return [p for p in paths if not os.path.exists(p)]


def boot_ci(v, n=2000, seed=0):
    rng = np.random.default_rng(seed)
    if len(v) == 0:
        return (np.nan, np.nan)
    bs = rng.choice(v, (n, len(v)), replace=True).mean(axis=1)
    return np.percentile(bs, 2.5), np.percentile(bs, 97.5)


def strip(ax, values, x, colour, rng, max_pts=140, width=0.09):
    v = values
    if len(v) > max_pts:
        v = rng.choice(v, max_pts, replace=False)
    ax.scatter(x + rng.normal(0, width, len(v)), v, s=3.5, alpha=0.28,
               color=colour, linewidths=0, zorder=2, rasterized=True)


# ------------------------------------------------------------------
def fig1_is_load(a):
    calls_p = os.path.join(a.isout, "IS_calls.tsv")
    miss = need(calls_p, a.meta)
    if miss:
        print(f"  skip fig1: missing {miss}")
        return
    calls = pd.read_csv(calls_p, sep="\t", low_memory=False)
    meta = pd.read_csv(a.meta, sep="\t", dtype={"fastbaps_group": str},
                       low_memory=False)
    meta["genome"] = meta.genome.astype(str)
    meta["fastbaps_group"] = meta.fastbaps_group.astype(str).str.strip()
    load = calls.groupby(calls.genome.astype(str)).size().rename("n_IS")
    d = meta.join(load, on="genome")
    d["n_IS"] = d.n_IS.fillna(0)
    if "total_len" in d.columns and d.total_len.notna().any():
        d["y"] = d.n_IS / d.total_len * 1e6
        ylab = "IS copies per Mb"
    else:
        d["y"] = d.n_IS
        ylab = "IS copies per genome"

    order = (d.groupby("fastbaps_group").y.median()
              .sort_values(ascending=False).index.tolist())
    has_forest = os.path.exists(os.path.join(a.isout,
                                             "02_focal_effect_summary.tsv"))
    fig, axes = plt.subplots(
        1, 2 if has_forest else 1,
        figsize=(7.2 if has_forest else 4.6, 3.1),
        gridspec_kw={"width_ratios": [1.75, 1]} if has_forest else None,
        constrained_layout=True)
    axes = np.atleast_1d(axes)
    ax = axes[0]

    rng = np.random.default_rng(0)
    for i, g in enumerate(order):
        v = d.loc[d.fastbaps_group == g, "y"].dropna().values
        strip(ax, v, i, FOCAL if g == a.focal else NEUTRAL, rng)
    bp = ax.boxplot([d.loc[d.fastbaps_group == g, "y"].dropna().values
                     for g in order],
                    positions=range(len(order)), widths=0.5,
                    patch_artist=True, showfliers=False, zorder=3,
                    medianprops=dict(color=INK, lw=1.3),
                    whiskerprops=dict(lw=0.7, color="#4D4D4D"),
                    capprops=dict(lw=0.7, color="#4D4D4D"))
    for patch, g in zip(bp["boxes"], order):
        focal_box = g == a.focal
        patch.set_facecolor(FOCAL_L if focal_box else "white")
        patch.set_edgecolor(FOCAL if focal_box else "#4D4D4D")
        patch.set_linewidth(1.0 if focal_box else 0.7)
        patch.set_alpha(0.95)

    ns = d.groupby("fastbaps_group").size()
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels(order)
    ax.set_xlabel("fastbaps group")
    ax.set_ylabel(ylab)
    ax.set_title(f"IS content across genetic groups", color=INK)
    top = ax.get_ylim()[1]
    ax.set_ylim(0, top * 1.12)
    for i, g in enumerate(order):
        ax.annotate(f"{ns[g]}", (i, top * 1.05), ha="center", fontsize=6,
                    color=FOCAL if g == a.focal else MUTED)
    ax.annotate("n", (-0.6, top * 1.05), ha="center", fontsize=6, color=MUTED)
    panel_letter(ax, "A")

    if has_forest:
        ax2 = axes[1]
        f = pd.read_csv(os.path.join(a.isout, "02_focal_effect_summary.tsv"),
                        sep="\t").iloc[::-1]
        y = np.arange(len(f))
        ax2.axvline(1, color=MUTED, lw=0.8, ls=(0, (3, 3)), zorder=1)
        ax2.hlines(y, f.ci_lo, f.ci_hi, color="#4D4D4D", lw=1.1, zorder=2)
        ax2.scatter(f.fold_change, y, s=26, color=FOCAL, zorder=3,
                    edgecolor="white", linewidth=0.6)
        ax2.set_yticks(y)
        ax2.set_yticklabels([m.split(":")[0].strip() for m in f.model])
        ax2.set_xlabel("fold change, focal vs other")
        ax2.set_title("effect across models", color=INK)
        ax2.set_ylim(-0.6, len(f) - 0.4)
        ax2.xaxis.set_major_locator(MaxNLocator(4))
        panel_letter(ax2, "B")

    save(fig, a.out, "fig1_IS_load")


# ------------------------------------------------------------------
def fig2_recency(a):
    """
    Recency of IS copies, built around the pair-count-adjusted estimate.

    The per-genome proportion p_recent is a proportion over n
    within-genome IS pairs, and n is small for most genomes (median ~4).
    Two consequences: the distribution spikes at simple fractions, and
    p_recent rises with n because more copies give more chances of a
    near-identical pair. The unweighted group difference therefore
    partly measures IS copy number, and in practice it changes sign as
    you filter on n. So:

      Panel A shows the distribution with genomes below --min-pairs
      excluded, on a proportion-of-genomes axis.
      Panel B shows the pair-weighted proportion, every genome counted
      in proportion to the evidence it carries.
      Panel C shows the odds ratio across model specifications, which is
      the number to report.
    """
    p = os.path.join(a.isout, "IS_recency_per_genome.tsv")
    if need(p):
        print(f"  skip fig2: missing {p}")
        return
    r = pd.read_csv(p, sep="\t").dropna(subset=["p_recent"])
    if "qfocal" not in r.columns:
        print("  skip fig2: no focal column")
        return
    has_n = "n" in r.columns
    n_before = len(r)
    r_f = r[r["n"] >= a.min_pairs] if has_n else r
    n_drop = n_before - len(r_f)

    sel = r_f.loc[r_f.qfocal == "focal", "p_recent"].values
    oth = r_f.loc[r_f.qfocal == "other", "p_recent"].values
    if len(sel) < 3 or len(oth) < 3:
        print(f"  skip fig2: too few genomes with >= {a.min_pairs} pairs")
        return

    eff_p = os.path.join(a.isout, "02_recency_effect_summary.tsv")
    has_eff = os.path.exists(eff_p)
    ncol = 3 if has_eff else 2
    fig, axes = plt.subplots(1, ncol, figsize=(3.3 * ncol, 3.1),
                             gridspec_kw={"width_ratios":
                                          [1.5, 0.9, 1.3][:ncol]},
                             constrained_layout=True)
    axes = np.atleast_1d(axes)

    ax = axes[0]
    bins = np.linspace(0, 1, 26)
    for v, c, lab in ((oth, NEUTRAL, f"other groups (n={len(oth)})"),
                      (sel, FOCAL, f"group {a.focal} (n={len(sel)})")):
        w = np.ones(len(v)) / len(v)
        ax.hist(v, bins=bins, weights=w, histtype="stepfilled", color=c,
                alpha=0.28, lw=0)
        ax.hist(v, bins=bins, weights=w, histtype="step", color=c, lw=1.5,
                label=lab)
    ax.set_xlabel("within-genome IS pairs >99% identical")
    ax.set_ylabel("proportion of genomes in group")
    ax.set_title("Recency of IS copies", color=INK)
    ax.legend(loc="upper left")
    ax.set_xlim(0, 1)
    if n_drop:
        ax.annotate(f"{n_drop} of {n_before} genomes excluded "
                    f"(<{a.min_pairs} pairs)", xy=(0, -0.21),
                    xycoords="axes fraction", fontsize=6.5, color=MUTED,
                    ha="left", va="top", annotation_clip=False)
    panel_letter(ax, "A")

    # pair-weighted proportion with a binomial CI: pairs, not genomes
    ax2 = axes[1]
    if has_n:
        for i, arm in enumerate(("other", "focal")):
            sub = r[r.qfocal == arm]
            k = (sub.p_recent * sub["n"]).round().sum()
            N = sub["n"].sum()
            if N == 0:
                continue
            ph = k / N
            se = np.sqrt(ph * (1 - ph) / N)
            c = NEUTRAL if arm == "other" else FOCAL
            ax2.vlines(i, ph - 1.96 * se, ph + 1.96 * se, color=INK, lw=1.5,
                       zorder=4)
            ax2.scatter(i, ph, s=40, color=c, zorder=5, edgecolor="white",
                        linewidth=0.8)
            ax2.annotate(f"{ph:.2f}\n({int(N):,} pairs)", (i, ph),
                         xytext=(9, 0), textcoords="offset points",
                         fontsize=6.5, color=MUTED, va="center")
        ax2.set_xticks([0, 1])
        ax2.set_xticklabels(["other", f"group {a.focal}"])
        ax2.set_xlim(-0.5, 1.8)
        ax2.set_ylabel("proportion of pairs >99% identical")
        ax2.set_title("Pair-weighted", color=INK)
    else:
        ax2.axis("off")
    panel_letter(ax2, "B")

    if has_eff:
        e = pd.read_csv(eff_p, sep="\t").iloc[::-1]
        ax3 = axes[2]
        y = np.arange(len(e))
        ax3.axvline(1, color=MUTED, lw=0.8, ls=(0, (3, 3)), zorder=1)
        ax3.hlines(y, e.ci_lo, e.ci_hi, color="#4D4D4D", lw=1.1, zorder=2)
        ax3.scatter(e.odds_ratio, y, s=26, color=FOCAL, zorder=3,
                    edgecolor="white", linewidth=0.6)
        ax3.set_yticks(y)
        ax3.set_yticklabels([m.split(":")[0].strip() for m in e.model],
                            fontsize=7.5)
        ax3.set_xlabel("odds ratio, focal vs other")
        ax3.set_title("Adjusted effect", color=INK)
        ax3.set_ylim(-0.6, len(e) - 0.4)
        ax3.xaxis.set_major_locator(MaxNLocator(4))
        panel_letter(ax3, "C")

    save(fig, a.out, "fig2_recency")

# ------------------------------------------------------------------
def fig3_sfs(a):
    p = os.path.join(a.isout, "insertion_frequencies.tsv")
    if need(p):
        print(f"  skip fig3: missing {p} (run script 03)")
        return
    f = pd.read_csv(p, sep="\t", index_col=0)
    if f.empty:
        print("  skip fig3: empty frequency table")
        return

    fig, axes = plt.subplots(1, 2, figsize=(7.0, 3.0), constrained_layout=True)
    ax = axes[0]
    ax.hist(f.freq, bins=np.linspace(0, 1, 26), color=FOCAL, alpha=0.85,
            edgecolor="white", linewidth=0.6)
    ax.set_xlabel("insertion frequency")
    ax.set_ylabel("insertion sites")
    ax.set_title(f"Site frequency spectrum (n = {len(f)})", color=INK)
    panel_letter(ax, "A")

    # the histogram hides the high-frequency tail when singletons dominate,
    # so show the cumulative form alongside it
    ax2 = axes[1]
    x = np.sort(f.freq.values)
    ax2.plot(x, np.arange(1, len(x) + 1) / len(x), color=FOCAL, lw=1.8)
    hi = float((f.freq >= 0.5).mean())
    ax2.axvline(0.5, color=MUTED, lw=0.8, ls=(0, (3, 3)))
    ax2.annotate(f"{hi:.1%} of sites\nat frequency ≥ 0.5",
                 xy=(0.5, 1 - hi), xytext=(0.12, max(0.18, 0.55 - hi)),
                 fontsize=7.5, color=INK,
                 arrowprops=dict(arrowstyle="->", lw=0.8, color=MUTED,
                                 shrinkB=3))
    ax2.set_xlabel("insertion frequency")
    ax2.set_ylabel("cumulative fraction of sites")
    ax2.set_title("Cumulative distribution", color=INK)
    ax2.set_xlim(0, 1)
    ax2.set_ylim(0, 1.02)
    panel_letter(ax2, "B")
    save(fig, a.out, "fig3_sfs")

    # ---- panel C: focal vs other groups, on a relative scale ----
    # Absolute counts cannot be compared: the focal group has more sites
    # and a different number of genomes. Recomputing frequency WITHIN each
    # group and plotting the proportion of that group's sites per bin puts
    # both on the same scale, which is what makes the shapes comparable.
    mat_p = os.path.join(a.isout, "insertion_matrix.tsv")
    if need(mat_p, a.meta):
        print("  fig3c skipped: needs insertion_matrix.tsv and metadata")
        return
    mat = pd.read_csv(mat_p, sep="\t", index_col=0, low_memory=False)
    meta = pd.read_csv(a.meta, sep="\t", dtype={"fastbaps_group": str},
                       low_memory=False)
    meta["genome"] = meta.genome.astype(str)
    meta["fastbaps_group"] = meta.fastbaps_group.astype(str).str.strip()
    gmap = dict(zip(meta.genome, meta.fastbaps_group))

    cols_focal = [c for c in mat.columns if gmap.get(str(c)) == a.focal]
    cols_other = [c for c in mat.columns
                  if str(c) in gmap and gmap[str(c)] != a.focal]
    if len(cols_focal) < 5 or len(cols_other) < 5:
        print("  fig3c skipped: too few genomes matched to groups")
        return

    def spectrum(cols, min_frac=0.5):
        sub = mat[cols]
        scored = sub.notna().sum(axis=1)
        keep = scored >= min_frac * len(cols)
        f = sub[keep].mean(axis=1, skipna=True).dropna()
        f = f[f > 0]                       # a site absent from the group is not segregating
        return f

    f_focal = spectrum(cols_focal)
    f_other = spectrum(cols_other)
    if len(f_focal) < 20 or len(f_other) < 20:
        print("  fig3c skipped: too few scorable sites per group")
        return

    bins = np.linspace(0, 1, 21)
    fig2, (axA, axB) = plt.subplots(1, 2, figsize=(7.2, 3.0),
                                    constrained_layout=True)
    w = (bins[1] - bins[0])
    centres = (bins[:-1] + bins[1:]) / 2
    hF, _ = np.histogram(f_focal, bins=bins)
    hO, _ = np.histogram(f_other, bins=bins)
    pF, pO = hF / hF.sum(), hO / hO.sum()
    axA.bar(centres - w / 4, pO, width=w / 2, color=NEUTRAL, alpha=0.85,
            label=f"other groups ({len(f_other)} sites, {len(cols_other)} genomes)")
    axA.bar(centres + w / 4, pF, width=w / 2, color=FOCAL, alpha=0.85,
            label=f"group {a.focal} ({len(f_focal)} sites, {len(cols_focal)} genomes)")
    axA.set_xlabel("insertion frequency within group")
    axA.set_ylabel("proportion of that group's sites")
    axA.set_title("Relative site frequency spectrum", color=INK)
    axA.legend(loc="upper right")
    panel_letter(axA, "A")

    # ratio view: where in the spectrum the groups differ
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = np.where(pO > 0, pF / pO, np.nan)
    axB.axhline(1, color=MUTED, lw=0.8, ls=(0, (3, 3)))
    ok = np.isfinite(ratio)
    axB.plot(centres[ok], ratio[ok], color=FOCAL, lw=1.6, marker="o",
             markersize=3.5, markeredgecolor="white", markeredgewidth=0.5)
    axB.set_yscale("log")
    axB.set_xlabel("insertion frequency within group")
    axB.set_ylabel(f"group {a.focal} / other groups")
    axB.set_title("Enrichment across the spectrum", color=INK)
    axB.annotate("above 1: relatively more\nsites at this frequency",
                 xy=(0.03, 0.95), xycoords="axes fraction", fontsize=6.5,
                 color=MUTED, va="top")
    panel_letter(axB, "B")

    save(fig2, a.out, "fig3c_relative_sfs")


# ------------------------------------------------------------------
def fig4_recombination(a):
    dens_p = os.path.join(a.isout, "recombination_density_by_group.tsv")
    corr_p = os.path.join(a.isout, "IS_vs_recombination_by_group.tsv")
    have = [p for p in (dens_p, corr_p) if os.path.exists(p)]
    if not have:
        print("  skip fig4: no recombination outputs (run 04b / 04c)")
        return

    n = len(have)
    fig, axes = plt.subplots(1, n, figsize=(3.7 * n, 3.0),
                             constrained_layout=True)
    axes = np.atleast_1d(axes)
    i = 0
    letters = iter("AB")

    if os.path.exists(dens_p):
        P = pd.read_csv(dens_p, sep="\t", index_col=0)
        ax = axes[i]; i += 1
        others = [c for c in P.columns if c != a.focal]
        if others:
            m = P[others].mean(axis=1)
            ax.fill_between(P.index / 1e6, 0, m, color=NEUTRAL_L, lw=0,
                            label="other groups (mean)")
            ax.plot(P.index / 1e6, m, color=NEUTRAL, lw=0.9)
        if a.focal in P.columns:
            ax.plot(P.index / 1e6, P[a.focal], color=FOCAL, lw=1.2,
                    label=f"group {a.focal}")
        ax.set_xlabel("core-genome position (Mb)")
        ax.set_ylabel("recombination density per genome")
        ax.set_title("Recombination along the core", color=INK)
        ax.legend(loc="upper right")
        ax.set_ylim(bottom=0)
        panel_letter(ax, next(letters))

    if os.path.exists(corr_p):
        c = pd.read_csv(corr_p, sep="\t").sort_values("spearman_rho")
        ax = axes[i]
        y = np.arange(len(c))
        ax.axvline(0, color=MUTED, lw=0.8, ls=(0, (3, 3)), zorder=1)
        cols = [FOCAL if str(g) == str(a.focal) else NEUTRAL for g in c.group]
        ax.hlines(y, 0, c.spearman_rho, color=cols, lw=1.0, alpha=0.5,
                  zorder=2)
        ax.scatter(c.spearman_rho, y, color=cols, s=26, zorder=3,
                   edgecolor="white", linewidth=0.6)
        ax.set_yticks(y)
        ax.set_yticklabels([f"{g}" for g in c.group], fontsize=7)
        ax.set_ylabel("fastbaps group")
        ax.set_xlabel("Spearman ρ, IS load vs recombination")
        ax.set_title("Within-group association", color=INK)
        ax.set_ylim(-0.7, len(c) - 0.3)
        panel_letter(ax, next(letters))

    save(fig, a.out, "fig4_recombination")


# ------------------------------------------------------------------
def fig5_selection(a):
    p = os.path.join(a.isout, "gene_distance_to_IS.tsv")
    if need(p):
        print(f"  skip fig5: missing {p} (run script 04)")
        return
    raw = pd.read_csv(p, sep="\t", low_memory=False)
    if not {"selected", "dist_IS"}.issubset(raw.columns):
        print("  skip fig5: unexpected columns")
        return
    raw["selected"] = raw.selected.astype(bool)
    raw = raw.replace([np.inf, -np.inf], np.nan).dropna(subset=["dist_IS"])

    # The table has one row per gene PER GENOME, so plotting it directly
    # counts each cluster once per genome carrying it: n is inflated by
    # roughly the number of genomes and the curves look far more precise
    # than the data warrant. Collapse to one value per cluster.
    if "cluster" in raw.columns:
        d = (raw.groupby("cluster")
                .agg(dist_IS=("dist_IS", "median"),
                     selected=("selected", "max"),
                     n_genomes=("dist_IS", "size"))
                .reset_index())
        unit = "accessory gene clusters"
        xlab = "median distance to nearest IS element (kb)"
    else:
        d = raw
        unit = "accessory genes"
        xlab = "distance to nearest IS element (kb)"
    d["selected"] = d.selected.astype(bool)

    fig, ax = plt.subplots(figsize=(4.6, 3.0), constrained_layout=True)
    for lab, c, name in ((False, NEUTRAL, "not selected"),
                         (True, FOCAL, "under selection")):
        v = d.loc[d.selected == lab, "dist_IS"].values
        if len(v) == 0:
            continue
        x = np.sort(v)
        ax.plot(x / 1000, np.arange(1, len(x) + 1) / len(x), color=c, lw=1.8,
                label=f"{name} (n={len(v)} genes)")
        ax.axvline(np.median(x) / 1000, color=c, lw=0.8, ls=(0, (2, 2)),
                   alpha=0.7)
    # symlog spans negatives by default, which wastes half the panel and
    # renders 10^0/10^1 mixed with a linear region; clamp and label plainly
    ax.set_xscale("symlog", linthresh=1, linscale=0.4)
    xmax = float(d.dist_IS.max()) / 1000
    ax.set_xlim(0, xmax * 1.15)
    ticks = [t for t in (0, 1, 10, 100, 1000) if t <= xmax * 1.15]
    ax.set_xticks(ticks)
    ax.set_xticklabels([str(t) for t in ticks])
    ax.set_xlabel(xlab)
    ax.set_ylabel(f"cumulative fraction of {unit}")
    ax.set_title("Proximity of accessory genes to IS elements", color=INK)
    ax.legend(loc="lower right")
    ax.set_ylim(0, 1.02)
    direction_note(ax, "← closer to IS", "further from IS →")
    save(fig, a.out, "fig5_selection")


# ------------------------------------------------------------------
def fig6_shell_is(a):
    """
    Shell-genome IS association, drawn to read the same way as fig5.

    Two deliberate choices so the reader never has to reverse their
    interpretation between figures:

      Panel A repeats fig5's exact construction (cumulative fraction of
      genes whose nearest IS lies within x kb) restricted to shell genes.
      Same axes, same statistic, different gene set.

      Panel B keeps the distinct quantity, the fraction of carrying
      genomes with an IS nearby, but as a COMPLEMENTARY cumulative: y is
      the fraction of genes IS-associated in at least x of their genomes.
      Drawn this way the more associated set is the HIGHER curve, as in
      panel A and fig5. The plain ECDF put it lower, which is what made
      the two figures look contradictory.
    """
    p = os.path.join(a.isout, "shell_cluster_IS_association.tsv")
    if need(p):
        print(f"  skip fig6: missing {p} (run 04g)")
        return
    d = pd.read_csv(p, sep="\t")
    if "prox_frac" not in d.columns or d.empty:
        print("  skip fig6: unexpected or empty table")
        return
    d["selected"] = d.selected.astype(bool)
    if d.selected.sum() < 3 or (~d.selected).sum() < 3:
        print("  skip fig6: too few clusters in one category")
        return
    shell_ids = set(d.cluster.astype(str)) if "cluster" in d.columns else set()

    dist_p = os.path.join(a.isout, "gene_distance_to_IS.tsv")
    have_dist = os.path.exists(dist_p) and bool(shell_ids)
    ncol = 3 if have_dist else 2
    fig, axes = plt.subplots(1, ncol, figsize=(3.4 * ncol, 3.1),
                             constrained_layout=True)
    axes = np.atleast_1d(axes)
    i = 0
    letters = iter("ABC")

    if have_dist:
        raw = pd.read_csv(dist_p, sep="\t", low_memory=False)
        raw = raw.replace([np.inf, -np.inf], np.nan).dropna(subset=["dist_IS"])
        raw["cluster"] = raw.cluster.astype(str)
        raw = raw[raw.cluster.isin(shell_ids)]
        if len(raw):
            g = (raw.groupby("cluster")
                    .agg(dist_IS=("dist_IS", "median"),
                         selected=("selected", "max")).reset_index())
            g["selected"] = g.selected.astype(bool)
            ax = axes[i]; i += 1
            for lab, c, name in ((False, NEUTRAL, "other shell genes"),
                                 (True, FOCAL, "under selection")):
                v = np.sort(g.loc[g.selected == lab, "dist_IS"].values)
                if len(v) == 0:
                    continue
                ax.plot(v / 1000, np.arange(1, len(v) + 1) / len(v), color=c,
                        lw=1.8, label=f"{name} (n={len(v)})")
            ax.set_xscale("symlog", linthresh=1, linscale=0.4)
            xmax = float(g.dist_IS.max()) / 1000
            ax.set_xlim(0, xmax * 1.15)
            ticks = [t for t in (0, 1, 10, 100, 1000) if t <= xmax * 1.15]
            ax.set_xticks(ticks)
            ax.set_xticklabels([str(t) for t in ticks])
            ax.set_xlabel("median distance to nearest IS (kb)")
            ax.set_ylabel("cumulative fraction of shell genes")
            ax.set_title("Distance to nearest IS", color=INK)
            ax.legend(loc="lower right")
            ax.set_ylim(0, 1.02)
            direction_note(ax, "\u2190 closer to IS", "further from IS \u2192")
            panel_letter(ax, next(letters))

    ax2 = axes[i]; i += 1
    for sub, c, lab in ((d[~d.selected], NEUTRAL, "other shell genes"),
                        (d[d.selected], FOCAL, "under selection")):
        x = np.sort(sub.prox_frac.values)
        y = 1 - np.arange(0, len(x)) / len(x)
        ax2.step(x, y, where="post", color=c, lw=1.8,
                 label=f"{lab} (n={len(sub)})")
    ax2.set_xlabel("threshold: fraction of carrying genomes")
    ax2.set_ylabel("fraction of shell genes at or above")
    ax2.set_title("IS-associated in at least x of genomes", color=INK)
    ax2.legend(loc="upper right")
    ax2.set_xlim(0, 1)
    ax2.set_ylim(0, 1.02)
    direction_note(ax2, "any association", "association in most genomes \u2192")
    panel_letter(ax2, next(letters))

    ax3 = axes[i]
    rng = np.random.default_rng(2)
    for j, (sub, c) in enumerate(((d[~d.selected], NEUTRAL),
                                  (d[d.selected], FOCAL))):
        v = sub.prox_frac.values
        strip(ax3, v, j, c, rng, max_pts=180, width=0.075)
        lo, hi = boot_ci(v)
        ax3.vlines(j, lo, hi, color=INK, lw=1.5, zorder=4)
        ax3.scatter(j, v.mean(), s=36, color=c, zorder=5,
                    edgecolor="white", linewidth=0.8)
    ax3.set_xticks([0, 1])
    ax3.set_xticklabels(["other", "selected"])
    ax3.set_xlim(-0.55, 1.55)
    ax3.set_ylabel("proximity fraction")
    ax3.set_title("Mean (95% CI)", color=INK)
    panel_letter(ax3, next(letters))

    save(fig, a.out, "fig6_shell_IS_frequency")


# ------------------------------------------------------------------
def fig7_neutral_sfs(a):
    """Insertion spectrum against the synonymous-SNP neutral reference."""
    p = os.path.join(a.isout, "sfs_projected.tsv")
    if need(p):
        print(f"  skip fig7: missing {p} (run 07_neutral_sfs_comparison.py)")
        return
    d = pd.read_csv(p, sep="\t")
    if not {"minor_allele_count", "synonymous_prop",
            "insertions_prop"}.issubset(d.columns):
        print("  skip fig7: unexpected columns")
        return

    cmp_p = os.path.join(a.isout, "sfs_masking_comparison.tsv")
    has_cmp = os.path.exists(cmp_p)
    n = 3 if has_cmp else 2
    fig, axes = plt.subplots(1, n, figsize=(3.5 * n, 3.1),
                             gridspec_kw={"width_ratios":
                                          [1.4, 1.2, 0.8][:n]},
                             constrained_layout=True)
    axes = np.atleast_1d(axes)

    k = d.minor_allele_count.values
    w = 0.42
    ax = axes[0]
    ax.bar(k - w / 2, d.synonymous_prop, width=w, color=NEUTRAL, alpha=0.9,
           label="synonymous SNPs")
    ax.bar(k + w / 2, d.insertions_prop, width=w, color=FOCAL, alpha=0.9,
           label="IS insertions")
    ax.set_yscale("log")
    ax.set_xlabel("minor allele count (projected)")
    ax.set_ylabel("proportion of sites")
    ax.set_title("Folded site frequency spectra", color=INK)
    ax.legend(loc="upper right")
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    panel_letter(ax, "A")

    # ratio: where in the spectrum insertions depart from neutrality
    ax2 = axes[1]
    with np.errstate(divide="ignore", invalid="ignore"):
        r = np.where(d.synonymous_prop > 0,
                     d.insertions_prop / d.synonymous_prop, np.nan)
    ok = np.isfinite(r) & (r > 0)
    ax2.axhline(1, color=MUTED, lw=0.8, ls=(0, (3, 3)))
    ax2.fill_between(k[ok], 1, r[ok], where=r[ok] > 1, color=FOCAL_L,
                     alpha=0.5, lw=0)
    ax2.plot(k[ok], r[ok], color=FOCAL, lw=1.6, marker="o", markersize=3.5,
             markeredgecolor="white", markeredgewidth=0.5)
    ax2.set_yscale("log")
    ax2.set_xlabel("minor allele count (projected)")
    ax2.set_ylabel("insertions / synonymous")
    ax2.set_title("Departure from neutrality", color=INK)
    ax2.xaxis.set_major_locator(MaxNLocator(integer=True))
    ax2.annotate("above 1: excess of insertions", xy=(0.03, 0.95),
                 xycoords="axes fraction", fontsize=6.5, color=MUTED,
                 va="top")
    panel_letter(ax2, "B")

    # masked vs unmasked reference: shows the result does not hinge on it
    if has_cmp:
        c = pd.read_csv(cmp_p, sep="\t")
        ax3 = axes[2]
        y = np.arange(len(c))
        ax3.axvline(1, color=MUTED, lw=0.8, ls=(0, (3, 3)))
        ax3.hlines(y, 1, c.enrichment, color=NEUTRAL, lw=1.0, alpha=0.6)
        ax3.scatter(c.enrichment, y, s=32, color=FOCAL, zorder=3,
                    edgecolor="white", linewidth=0.6)
        ax3.set_yticks(y)
        ax3.set_yticklabels(c.reference)
        ax3.set_xlabel("high-frequency enrichment")
        ax3.set_title("Neutral reference", color=INK)
        ax3.set_ylim(-0.6, len(c) - 0.4)
        for yi, (e, ns) in enumerate(zip(c.enrichment, c.n_syn_sites)):
            ax3.annotate(f"{e:.1f}x  ({ns} sites)", (e, yi), xytext=(5, 0),
                         textcoords="offset points", fontsize=6.5,
                         va="center", color=MUTED)
        ax3.set_xlim(0, c.enrichment.max() * 1.6)
        panel_letter(ax3, "C")

    save(fig, a.out, "fig7_neutral_sfs")


FIGS = {"fig1": fig1_is_load, "fig2": fig2_recency, "fig3": fig3_sfs,
        "fig4": fig4_recombination, "fig5": fig5_selection,
        "fig6": fig6_shell_is, "fig7": fig7_neutral_sfs}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--isout", default="is_out")
    ap.add_argument("--meta", default="metadata.tsv")
    ap.add_argument("--focal", required=True)
    ap.add_argument("--out", default="figures")
    ap.add_argument("--min-pairs", type=int, default=10,
                    help="fig2: minimum within-genome IS pairs for a genome "
                         "to contribute (0 keeps everything)")
    ap.add_argument("--only", nargs="*", default=None,
                    help="e.g. --only fig3 fig6")
    a = ap.parse_args()
    a.focal = str(a.focal).strip()
    os.makedirs(a.out, exist_ok=True)
    print(f"writing figures to {a.out}/")
    for key, fn in FIGS.items():
        if a.only and key not in a.only:
            continue
        try:
            fn(a)
        except Exception as e:
            print(f"  {key} failed: {type(e).__name__}: {e}")
    print("done")


if __name__ == "__main__":
    main()
