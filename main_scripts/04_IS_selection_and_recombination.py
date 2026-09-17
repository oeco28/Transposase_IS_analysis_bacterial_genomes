#!/usr/bin/env python3
"""
04_IS_selection_and_recombination.py

Three tests that connect the IS elements to the biology you care about.

  A. Physical association: are accessory genes under selection closer to
     IS elements, and more often flanked on both sides (composite
     transposon architecture), than accessory genes that are not under
     selection? Null must control for accessory-genome frequency, because
     rare genes are both more likely to be recent arrivals and more
     likely to sit in MGE-rich regions.

  B. NOTE: test_B below assumes ONE gubbins run against ONE reference.
     If gubbins was run per fastbaps group on a CORE alignment (the
     usual recommendation), use 04b_gubbins_integration.py instead:
     IS elements are mostly in the accessory genome, which a core
     alignment excludes, and per-group runs have no shared coordinate
     system. test_B remains valid only for the minority of IS copies
     inside core regions, run separately per group, never pooled.

  B. Transposition vs homologous recombination: two complementary
     questions, if you do have whole-genome single-reference blocks.
       B1. Do IS insertions fall inside gubbins recombination blocks?
           If the IS-associated accessory gene gains are largely OUTSIDE
           recombination tracts, transposition is acting as an
           independent force rather than being cargo of recombination.
       B2. Are gubbins block BOUNDARIES enriched for IS elements?
           If so, IS copies are also templating homologous exchange
           between genomes (ectopic recombination), which is a second,
           distinct mechanism by which the same elements inflate HGT in
           this group. Either result is publishable; both together are
           a much stronger story.

  C. Parallelism: independent insertions at the same locus on different
     branches of the tree. Convergence is the single most persuasive
     signal that an insertion is adaptive rather than drifting, and it
     is the analysis reviewers will find hardest to argue with.

Inputs
  is_out/IS_calls.tsv                per-copy IS table (script 01)
  gff/<genome>.gff                   Prokka/Bakta annotation
  gene_presence_absence.csv          Panaroo/Roary
  selected_genes.txt                 one gene cluster name per line
  gubbins_recombination_predictions.gff
  core_gene_tree.nwk
"""

import argparse, random, re
from collections import defaultdict
import numpy as np, pandas as pd
from scipy import stats

# ------------------------------------------------------------------
# A. distance from selected accessory genes to nearest IS
# ------------------------------------------------------------------

def parse_gff(path):
    rows = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                if line.startswith("##FASTA"):
                    break
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "CDS":
                continue
            gid = re.search(r"ID=([^;]+)", f[8])
            rows.append((f[0], int(f[3]), int(f[4]), f[6],
                         gid.group(1) if gid else None))
    return pd.DataFrame(rows, columns=["contig", "start", "end", "strand", "locus_tag"])


def nearest_is(genes, is_tab):
    """Distance from each CDS to the nearest IS on the same contig."""
    out = []
    for contig, gsub in genes.groupby("contig"):
        isc = is_tab[is_tab.contig == contig]
        if isc.empty:
            out.append(pd.Series(np.inf, index=gsub.index))
            continue
        s, e = isc["start"].values, isc["end"].values
        d = []
        for _, g in gsub.iterrows():
            left = g["start"] - e
            right = s - g["end"]
            overlap = ((s <= g["end"]) & (e >= g["start"])).any()
            d.append(0 if overlap else
                     min(np.min(left[left > 0], initial=np.inf),
                         np.min(right[right > 0], initial=np.inf)))
        out.append(pd.Series(d, index=gsub.index))
    return pd.concat(out).reindex(genes.index)


def flanked_both_sides(genes, is_tab, window=15000):
    """Composite-transposon-like: IS within `window` on both sides."""
    flag = {}
    for contig, gsub in genes.groupby("contig"):
        isc = is_tab[is_tab.contig == contig]
        s, e = isc["start"].values, isc["end"].values
        for i, g in gsub.iterrows():
            up = ((e <= g["start"]) & (e > g["start"] - window)).any()
            dn = ((s >= g["end"]) & (s < g["end"] + window)).any()
            flag[i] = bool(up and dn)
    return pd.Series(flag).reindex(genes.index)


def test_A(is_calls, gff_dir, pa_path, selected_path, focal_genomes, nperm=10000):
    pa = pd.read_csv(pa_path, low_memory=False)
    locus_to_cluster = {}
    idcols = [c for c in pa.columns if c in focal_genomes]
    for _, r in pa.iterrows():
        for c in idcols:
            v = r[c]
            if isinstance(v, str):
                for lt in re.split(r"[;\t ]+", v):
                    if lt:
                        locus_to_cluster[lt] = r.iloc[0]
    selected = set(open(selected_path).read().split())
    # accessory-frequency strata for the null: compare like with like
    freq = pa[idcols].notna().sum(axis=1) / len(idcols)
    cluster_freq = dict(zip(pa.iloc[:, 0], freq))

    recs = []
    for g in focal_genomes:
        genes = parse_gff(f"{gff_dir}/{g}.gff")
        ist = is_calls[is_calls.genome == g]
        if genes.empty:
            continue
        genes["dist_IS"] = nearest_is(genes, ist).values
        genes["flanked"] = flanked_both_sides(genes, ist).values
        genes["cluster"] = genes.locus_tag.map(locus_to_cluster)
        genes["genome"] = g
        recs.append(genes)
    G = pd.concat(recs)
    G = G.dropna(subset=["cluster"])
    G["freq"] = G.cluster.map(cluster_freq)
    G["accessory"] = G.freq < 0.95
    G["selected"] = G.cluster.isin(selected)

    acc = G[G.accessory & np.isfinite(G.dist_IS)]
    obs = acc.loc[acc.selected, "dist_IS"].median()

    # stratified permutation: shuffle the "selected" label within
    # deciles of accessory-genome frequency, preserving how common the
    # genes are, so the result cannot be explained by rare genes simply
    # living in MGE-dense neighbourhoods
    acc = acc.copy()
    acc["bin"] = pd.qcut(acc.freq, 10, duplicates="drop")
    null = []
    for _ in range(nperm):
        lab = acc.groupby("bin", observed=True)["selected"].transform(
            lambda x: np.random.permutation(x.values))
        null.append(acc.loc[lab.astype(bool), "dist_IS"].median())
    p = (np.sum(np.array(null) <= obs) + 1) / (nperm + 1)
    print(f"[A] median distance to nearest IS, selected accessory genes: "
          f"{obs:.0f} bp; stratified permutation p = {p:.4g}")

    ct = pd.crosstab(acc.selected, acc.flanked)
    print("[A] composite-transposon flanking, Fisher:",
          stats.fisher_exact(ct))

    # persist for 05_figures.py (fig5) and for reporting exact numbers
    import os
    os.makedirs("is_out", exist_ok=True)
    acc[["genome", "cluster", "freq", "selected", "flanked", "dist_IS"]].to_csv(
        "is_out/gene_distance_to_IS.tsv", sep="\t", index=False)
    with open("is_out/04_testA_summary.tsv", "w") as fh:
        fh.write("statistic\tvalue\n")
        fh.write(f"median_dist_selected\t{obs}\n")
        fh.write(f"median_dist_not_selected\t"
                 f"{acc.loc[~acc.selected,'dist_IS'].median()}\n")
        fh.write(f"permutation_p\t{p}\n")
        fh.write(f"n_selected_genes\t{int(acc.selected.sum())}\n")
        fh.write(f"n_other_accessory\t{int((~acc.selected).sum())}\n")
        fr = stats.fisher_exact(ct)
        fh.write(f"flanking_odds_ratio\t{fr[0]}\n")
        fh.write(f"flanking_p\t{fr[1]}\n")
    print("[A] wrote is_out/gene_distance_to_IS.tsv and "
          "is_out/04_testA_summary.tsv")
    return G


# ------------------------------------------------------------------
# B. IS vs gubbins recombination blocks
# ------------------------------------------------------------------

def parse_gubbins(path):
    rows = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t")
            if len(f) < 9:
                continue
            taxa = re.search(r'taxa="([^"]+)"', f[8])
            rows.append((int(f[3]), int(f[4]),
                         taxa.group(1).split() if taxa else []))
    return pd.DataFrame(rows, columns=["start", "end", "taxa"])


def test_B(is_ref_coords, gubbins_path, ref_len, nperm=10000, win=2000):
    """is_ref_coords: IS positions projected onto the gubbins reference."""
    gb = parse_gubbins(gubbins_path)
    breakpoints = np.concatenate([gb.start.values, gb.end.values])

    pos = np.asarray(is_ref_coords)
    inside = np.array([((gb.start <= p) & (gb.end >= p)).any() for p in pos])
    frac_block = (gb.end - gb.start).sum() / ref_len
    print(f"[B1] IS inside recombination blocks: {inside.mean():.3f} "
          f"(genome fraction covered by blocks: {frac_block:.3f}); "
          f"binomial p = {stats.binomtest(inside.sum(), len(pos), frac_block).pvalue:.3g}")

    d = np.array([np.min(np.abs(breakpoints - p)) for p in pos])
    obs = np.mean(d < win)
    null = []
    for _ in range(nperm):
        rp = np.random.randint(0, ref_len, size=len(pos))
        dn = np.array([np.min(np.abs(breakpoints - p)) for p in rp])
        null.append(np.mean(dn < win))
    p = (np.sum(np.array(null) >= obs) + 1) / (nperm + 1)
    print(f"[B2] IS within {win} bp of a recombination breakpoint: "
          f"{obs:.3f} vs null {np.mean(null):.3f}, p = {p:.4g}")


# ------------------------------------------------------------------
# C. parallel (convergent) insertions
# ------------------------------------------------------------------

def test_C(matrix_path, tree_path, gene_of_site=None):
    """
    Count minimal independent gains per insertion site by parsimony on
    the core tree. A site gained >=2 times independently is a strong
    adaptive candidate; a GENE hit by independent insertions at
    different sites in different lineages is stronger still.
    Requires: pip install pastml  (or use ape/phangorn in R)
    """
    from ete3 import Tree
    mat = pd.read_csv(matrix_path, sep="\t", index_col=0)
    t = Tree(tree_path, format=1)
    results = {}
    for site, row in mat.iterrows():
        states = {k: int(v) for k, v in row.items() if not pd.isna(v)}
        if len(states) < 10 or sum(states.values()) == 0:
            continue
        gains = 0
        for node in t.traverse("postorder"):
            if node.is_leaf():
                node.st = states.get(node.name, None)
            else:
                ch = [c.st for c in node.children if c.st is not None]
                node.st = 1 if ch and all(c == 1 for c in ch) else (0 if ch else None)
        for node in t.traverse("preorder"):
            if node.st == 1 and (node.up is None or node.up.st == 0):
                gains += 1
        results[site] = gains
    r = pd.Series(results).sort_values(ascending=False)
    print("[C] sites with >=2 independent gains:", (r >= 2).sum())
    r.to_csv("is_out/independent_gains_per_site.tsv", sep="\t")
    return r


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--calls", default="is_out/IS_calls.tsv")
    ap.add_argument("--gff", default="gff")
    ap.add_argument("--pa", default="gene_presence_absence.csv")
    ap.add_argument("--selected", default="selected_genes.txt")
    ap.add_argument("--meta", default="metadata.tsv")
    ap.add_argument("--focal", default="G4")
    a = ap.parse_args()

    meta = pd.read_csv(a.meta, sep="\t")
    focal = meta.loc[meta.fastbaps_group == a.focal, "genome"].astype(str).tolist()
    is_calls = pd.read_csv(a.calls, sep="\t")
    test_A(is_calls, a.gff, a.pa, a.selected, focal)
