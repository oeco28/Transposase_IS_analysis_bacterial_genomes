#!/usr/bin/env python3
"""
group_loci.py
Is the IS-proximity signal for genes under selection (focal group) driven by a few loci?

Neighbouring selected genes are not independent observations: one island of 8 selected genes next
to one IS counts as 8 genes in 01_is_proximity.py. This script
  1. builds a PER-GENE table for the focal group (which selected genes are IS-proximal?),
  2. collapses genes into LOCI (accessory genes that sit within --locus_gap bp of each other
     in the same order in >= --min_frac of the focal-group genomes that carry both),
  3. repeats the permutation test with the LOCUS as the unit
     (label = locus contains >= 1 selected gene; value = mean fraction of copies within --window of an IS;
      labels permuted within strata of locus frequency),
  4. runs a leave-one-locus-out check (does the result depend on one island?).

IS-encoded gene clusters are excluded exactly as in 01_is_proximity.py.
Inputs are the same files as 01_is_proximity.py (cds, is_table, genes, groups).

Outputs (in --out)
  group<focal>_genes.tsv            per-gene: selected, prop_near, med_dist, locus_id
  group<focal>_loci.tsv             per-locus: genes, n_selected, mean prop_near, selected_locus
  group<focal>_locus_tests.tsv      gene-level vs locus-level test, for every --locus_gaps value
  group<focal>_leave_one_out.tsv    effect/p after dropping each selected locus (first gap only)
"""
import argparse
import importlib.util
import os
import sys

import numpy as np
import pandas as pd
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from utils import make_strata, stratified_permutation  # noqa: E402

_spec = importlib.util.spec_from_file_location("isprox", os.path.join(HERE, "01_is_proximity.py"))
isprox = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(isprox)


def build_loci(cds_focal, gene_n, gap, min_pair_genomes, min_frac):
    """Connected components of accessory genes that are adjacent (<= gap bp) in the focal group."""
    d = cds_focal.sort_values(["genome", "contig", "start"])
    g, c = d["genome"].to_numpy(), d["contig"].to_numpy()
    s, e = d["start"].to_numpy(), d["end"].to_numpy()
    gc = d["gene_cluster"].to_numpy()
    n, found = len(d), []
    # pair every accessory gene with ALL later accessory genes within `gap` bp on the same contig
    # (not only the next one), so interloper genes in some genomes do not break a locus
    for k in range(1, min(n, 500)):
        same = (g[:-k] == g[k:]) & (c[:-k] == c[k:]) & ((s[k:] - e[:-k]) <= gap)
        if not same.any():
            break
        ok = same & (gc[:-k] != gc[k:])
        a_, b_ = gc[:-k][ok], gc[k:][ok]
        found.append(pd.DataFrame({"a": np.minimum(a_, b_), "b": np.maximum(a_, b_), "genome": g[:-k][ok]}))
    pr = (pd.concat(found).drop_duplicates() if found
          else pd.DataFrame({"a": [], "b": [], "genome": []}))
    cnt = pr.groupby(["a", "b"]).size().rename("n_pair").reset_index()
    cnt["frac"] = cnt["n_pair"] / np.minimum(cnt["a"].map(gene_n), cnt["b"].map(gene_n))
    cnt = cnt[(cnt["n_pair"] >= min_pair_genomes) & (cnt["frac"] >= min_frac)]
    genes = pd.Index(sorted(gene_n.index))
    idx = {x: i for i, x in enumerate(genes)}
    n = len(genes)
    if len(cnt):
        ia, ib = cnt["a"].map(idx).to_numpy(), cnt["b"].map(idx).to_numpy()
        m = coo_matrix((np.ones(len(ia)), (ia, ib)), shape=(n, n))
    else:
        m = coo_matrix((n, n))
    _, lab = connected_components(m, directed=False)
    return pd.Series(["L%04d" % x for x in lab], index=genes, name="locus_id")


def locus_table(gtab):
    L = gtab.groupby("locus_id").agg(
        n_genes=("selected", "size"), n_selected=("selected", "sum"),
        mean_prop_near=("prop_near", "mean"), max_prop_near=("prop_near", "max"),
        mean_local_freq=("local_freq", "mean"),
        genes=("gene_cluster", lambda x: ";".join(x)),
        products=("product", lambda x: " | ".join(sorted(set(str(v) for v in x)))[:120]))
    L["selected_locus"] = (L["n_selected"] > 0).astype(int)
    return L.reset_index()


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cds", required=True)
    ap.add_argument("--is_table", required=True)
    ap.add_argument("--genes", required=True)
    ap.add_argument("--groups", required=True)
    ap.add_argument("--focal", required=True)
    ap.add_argument("--out", default="results_group_loci")
    ap.add_argument("--window", type=int, default=5000)
    ap.add_argument("--censor", type=int, default=50000)
    ap.add_argument("--nperm", type=int, default=10000)
    ap.add_argument("--loo_nperm", type=int, default=2000)
    ap.add_argument("--min_genomes", type=int, default=5)
    ap.add_argument("--locus_gaps", default="10000", help="comma list of bp; first one is used for detail files")
    ap.add_argument("--min_pair_genomes", type=int, default=3)
    ap.add_argument("--min_frac", type=float, default=0.5)
    ap.add_argument("--is_overlap_frac", type=float, default=0.5)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    focal = str(a.focal)

    cols = {"genome": str, "contig": str, "gene_cluster": str}
    cds = isprox.norm_coords(pd.read_csv(a.cds, sep="\t", dtype=cols))
    ist = isprox.norm_coords(pd.read_csv(a.is_table, sep="\t", dtype={"genome": str, "contig": str}))
    genes = pd.read_csv(a.genes, sep="\t", dtype={"gene_cluster": str})
    groups = pd.read_csv(a.groups, sep="\t", dtype={"genome": str, "group": str})
    if "complete" in ist.columns:
        ist = ist[ist["complete"].astype(str).str.lower().isin(["1", "true", "t", "yes"])]

    # focal-group genomes only
    fg = set(groups.loc[groups["group"] == focal, "genome"])
    cds = cds[cds["genome"].isin(fg)].copy()
    ist = ist[ist["genome"].isin(fg)]
    print(f"focal group {focal}: {len(fg)} genomes, {len(cds)} CDS, {len(ist)} complete IS")

    cds["dist"] = isprox.nearest_distance(cds, ist, a.censor)
    cds["overlap_is"] = cds["dist"] == 0
    cl = cds.groupby("gene_cluster")["overlap_is"].mean().rename("frac_overlap_is")
    genes = genes.merge(cl, left_on="gene_cluster", right_index=True, how="left")
    prod = genes["product"].fillna("") if "product" in genes.columns else pd.Series("", index=genes.index)
    genes["is_encoded"] = ((genes["frac_overlap_is"].fillna(0) >= a.is_overlap_frac)
                           | prod.str.contains(isprox.IS_PRODUCT_PATTERN, regex=True))
    keep = genes[(genes["class"].str.lower() != "core") & (~genes["is_encoded"])].set_index("gene_cluster")
    cds = cds[cds["gene_cluster"].isin(keep.index)].copy()
    cds["near"] = cds["dist"] <= a.window

    s = cds.groupby("gene_cluster").agg(n_genomes=("genome", "nunique"), prop_near=("near", "mean"),
                                        med_dist=("dist", "median"))
    s = s.join(keep[["selected"] + (["product"] if "product" in keep.columns else [])], how="inner")
    if "product" not in s.columns:
        s["product"] = ""
    s = s[s["n_genomes"] >= a.min_genomes].copy()
    s["local_freq"] = s["n_genomes"] / len(fg)
    s["selected"] = s["selected"].astype(int)
    cds = cds[cds["gene_cluster"].isin(s.index)]
    print(f"accessory genes tested: {len(s)}; selected: {int(s.selected.sum())}")

    sel = s[s.selected == 1]
    print("\nSELECTED GENES IN FOCAL GROUP:")
    print(f"  with any copy within {a.window} bp of an IS : {(sel.prop_near > 0).sum()} / {len(sel)}")
    print(f"  with >= 50% of copies within the window      : {(sel.prop_near >= 0.5).sum()} / {len(sel)}")
    print(f"  mean fraction of copies near IS (selected)    : {sel.prop_near.mean():.3f}; other: "
          f"{s[s.selected == 0].prop_near.mean():.3f}")

    # gene-level test (same as 01) for reference
    gene_res = stratified_permutation(s["prop_near"], s["selected"].astype(bool),
                                      make_strata(s["local_freq"]), nperm=a.nperm)

    rows = []
    gaps = [int(x) for x in a.locus_gaps.split(",")]
    for k, gap in enumerate(gaps):
        locus_of = build_loci(cds, s["n_genomes"], gap, a.min_pair_genomes, a.min_frac)
        g2 = s.join(locus_of, how="left")
        g2.index.name = "gene_cluster"
        g2 = g2.reset_index()
        g2["locus_id"] = g2["locus_id"].fillna("S_" + g2["gene_cluster"].astype(str))
        L = locus_table(g2)
        res = stratified_permutation(L["mean_prop_near"], L["selected_locus"].astype(bool),
                                     make_strata(L["mean_local_freq"]), nperm=a.nperm)
        n_sel_loci = int(L["selected_locus"].sum())
        n_multi = int(((L["selected_locus"] == 1) & (L["n_selected"] > 1)).sum())
        rows.append(dict(locus_gap=gap, n_loci=len(L), n_selected_genes=int(s.selected.sum()),
                         n_selected_loci=n_sel_loci, selected_loci_with_multiple_selected_genes=n_multi,
                         gene_level_diff=gene_res["obs"], gene_level_p=gene_res["p_greater"],
                         locus_level_diff=res["obs"], locus_level_p=res["p_greater"],
                         locus_prop_near_selected=res["mean_sel"], locus_prop_near_other=res["mean_other"]))
        if k == 0:
            g2.to_csv(f"{a.out}/group{focal}_genes.tsv", sep="\t", index=False)
            L.sort_values(["selected_locus", "mean_prop_near"], ascending=False).to_csv(
                f"{a.out}/group{focal}_loci.tsv", sep="\t", index=False)
            # leave-one-locus-out
            loo = []
            for lid in L.loc[L.selected_locus == 1, "locus_id"]:
                Lr = L[L.locus_id != lid]
                r = stratified_permutation(Lr["mean_prop_near"], Lr["selected_locus"].astype(bool),
                                           make_strata(Lr["mean_local_freq"]), nperm=a.loo_nperm)
                loo.append(dict(dropped_locus=lid,
                                n_selected_genes_in_locus=int(L.loc[L.locus_id == lid, "n_selected"].iloc[0]),
                                locus_mean_prop_near=float(L.loc[L.locus_id == lid, "mean_prop_near"].iloc[0]),
                                diff=r["obs"], p_greater=r["p_greater"]))
            loo = pd.DataFrame(loo).sort_values("p_greater", ascending=False)
            loo.to_csv(f"{a.out}/group{focal}_leave_one_out.tsv", sep="\t", index=False)
            print(f"\nLeave-one-locus-out (gap {gap} bp): diff range {loo['diff'].min():+.4f} to "
                  f"{loo['diff'].max():+.4f}; worst-case p = {loo['p_greater'].max():.4g} "
                  f"(dropping {loo.iloc[0]['dropped_locus']})")
    out = pd.DataFrame(rows)
    out.to_csv(f"{a.out}/group{focal}_locus_tests.tsv", sep="\t", index=False)
    print("\nGENE-LEVEL vs LOCUS-LEVEL permutation tests (one-sided):")
    print(out.round(4).to_string(index=False))


if __name__ == "__main__":
    main()
