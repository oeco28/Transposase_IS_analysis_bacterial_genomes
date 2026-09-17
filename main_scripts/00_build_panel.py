#!/usr/bin/env python3
"""
00_build_panel.py  (v2)

Decides WHICH genomes go into which stage. Run before script 01.

  DISCOVERY (script 01, ISEScan)          -> ALL genomes, no exceptions.
     Cheap, and restricting it creates ascertainment bias in the site
     frequency spectrum that you cannot undo afterwards.

  SCORING   (script 03, insertion matrix) -> ALL focal genomes plus a
     stratified panel of non-focal genomes. Cost scales as
     (IS copies x genomes), so this is the only stage where subsetting
     is defensible. You still need enough genomes in each sister group
     to observe the EMPTY allele that polarises an insertion.

v2 changes: works with any metadata file rather than only the one
produced by 00a_prep_inputs.sh. Column names are auto-detected or given
explicitly with --id-col / --group-col. Group labels are coerced to
strings so integer fastbaps labels (--focal 12) match. Tree/metadata ID
overlap is checked and reported instead of silently yielding an empty
frame.

Usage:
  python 00_build_panel.py --meta my_metadata.tsv --tree tree.nwk --focal 12
  python 00_build_panel.py --meta my_metadata.tsv --tree tree.nwk --focal 12 \
      --id-col isolate_id --group-col fastbaps_level2
  python 00_build_panel.py --meta my_metadata.tsv --list-columns
"""

import argparse
import re
import sys

import numpy as np
import pandas as pd
from ete3 import Tree

ID_PAT = r"^(genome|isolate|sample|strain|taxon|accession|assembly|id|name|" \
         r".*_id|isolate.*|genome.*|sample.*)$"
GROUP_PAT = r"(fastbaps|baps|cluster|group|partition|lineage|pop)"
# Columns that look like group columns but are not. "Assembly Level" from an
# NCBI Datasets table is the one that bites: it matched the old pattern.
GROUP_EXCLUDE = r"(assembly|annotation|checkm|organism|wgs|accession|date|tech)"


def sniff_sep(path):
    if path.endswith((".tsv", ".txt", ".tab")):
        return "\t"
    if path.endswith(".csv"):
        return ","
    with open(path) as fh:
        head = fh.readline()
    return "\t" if head.count("\t") > head.count(",") else ","


def pick_column(cols, pattern, kind, explicit=None):
    if explicit is not None:
        if explicit not in cols:
            sys.exit(f"--{kind}-col '{explicit}' is not in the file.\n"
                     f"Columns present: {list(cols)}")
        return explicit
    hits = [c for c in cols if re.search(pattern, str(c), re.I)]
    if kind == "group":
        hits = [c for c in hits if not re.search(GROUP_EXCLUDE, str(c), re.I)]
    if not hits:
        sys.exit(f"Could not find a {kind} column automatically.\n"
                 f"Columns present: {list(cols)}\n"
                 f"Rerun with --{kind}-col <name>")
    if len(hits) > 1:
        print(f"note: several candidate {kind} columns {hits}; using '{hits[0]}'. "
              f"Override with --{kind}-col if that is wrong.")
    return hits[0]


def strip_ext(s):
    return (s.astype(str)
             .str.strip()
             .str.replace(r"\.(fna|fa|fasta|gbk|gff)(\.gz)?$", "", regex=True))


def dereplicate(tree, thresh):
    """
    Maximal clades whose maximum internal pairwise distance is below
    `thresh`, treated as clonal clusters. Bounded by 2 x the deepest
    root-to-leaf distance within the clade, which is cheap and slightly
    conservative (it may split a cluster, never merge two).
    """
    for n in tree.traverse("postorder"):
        n.maxdepth = 0.0 if n.is_leaf() else max(c.maxdepth + c.dist
                                                 for c in n.children)
    clusters, stack = [], [tree]
    while stack:
        node = stack.pop()
        if node.is_leaf():
            clusters.append([node.name])
        elif 2 * node.maxdepth <= thresh:
            clusters.append(node.get_leaf_names())
        else:
            stack.extend(node.children)
    return clusters


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--meta", required=True)
    ap.add_argument("--tree", default="core_gene_tree.nwk")
    ap.add_argument("--focal", required=True,
                    help="focal fastbaps group label, e.g. 12")
    ap.add_argument("--id-col", default=None)
    ap.add_argument("--group-col", default=None)
    ap.add_argument("--list-columns", action="store_true",
                    help="print the metadata columns and exit")
    ap.add_argument("--derep-thresh", type=float, default=1e-4,
                    help="patristic distance below which tips are called clonal")
    ap.add_argument("--max-nonfocal", type=int, default=300,
                    help="cap on non-focal genomes in the SCORING panel only")
    ap.add_argument("--out-prefix", default="panel")
    a = ap.parse_args()

    meta = pd.read_csv(a.meta, sep=sniff_sep(a.meta), dtype=str,
                       low_memory=False)
    if a.list_columns:
        print("\n".join(f"  {c}" for c in meta.columns))
        return

    id_col = pick_column(meta.columns, ID_PAT, "id", a.id_col)
    grp_col = pick_column(meta.columns, GROUP_PAT, "group", a.group_col)
    print(f"genome ID column : {id_col}")
    print(f"group column     : {grp_col}")

    meta = meta.rename(columns={id_col: "genome", grp_col: "fastbaps_group"})
    meta["genome"] = strip_ext(meta["genome"])
    # coerce to string on both sides so --focal 12 matches an integer column
    meta["fastbaps_group"] = meta["fastbaps_group"].astype(str).str.strip()
    focal_label = str(a.focal).strip()

    groups = meta.fastbaps_group.value_counts()
    if focal_label not in set(groups.index):
        sys.exit(f"\n--focal '{focal_label}' is not a value in '{grp_col}'.\n"
                 f"Values present: {list(groups.index)[:30]}")
    print(f"\n{len(meta)} genomes in metadata, "
          f"{meta.fastbaps_group.nunique()} groups; "
          f"focal group '{focal_label}' has {groups[focal_label]}")

    # ---- tree / metadata ID reconciliation ----------------------
    tree = Tree(a.tree, format=1)
    tips = pd.Series(tree.get_leaf_names())
    tips_clean = set(strip_ext(tips))
    overlap = set(meta.genome) & tips_clean
    if len(overlap) == 0:
        sys.exit("\nNo genome IDs are shared between the metadata and the tree.\n"
                 f"  metadata examples: {meta.genome.head(3).tolist()}\n"
                 f"  tree tip examples: {list(tips_clean)[:3]}\n"
                 "These are formatted differently. Fix the IDs, or pass the "
                 "correct column with --id-col.")
    if len(overlap) < 0.9 * min(len(meta), len(tips_clean)):
        print(f"\nWARNING: only {len(overlap)} IDs shared "
              f"({len(meta)} in metadata, {len(tips_clean)} tips). "
              "Unmatched genomes are dropped from every downstream model, "
              "so check this is intended.")
        missing = sorted(set(meta.genome) - tips_clean)[:5]
        if missing:
            print(f"  in metadata but not the tree, e.g. {missing}")

    meta = meta[meta.genome.isin(tips_clean)].copy()
    tree.prune([t for t in tree.get_leaf_names()
                if re.sub(r"\.(fna|fa|fasta)(\.gz)?$", "", t) in set(meta.genome)],
               preserve_branch_length=True)
    print(f"{len(meta)} genomes retained after matching to the tree")

    # ---- discovery panel: everything ----------------------------
    meta.genome.to_csv(f"{a.out_prefix}_discovery.txt", index=False, header=False)
    print(f"\ndiscovery panel: {len(meta)} genomes (all of them) -> "
          f"{a.out_prefix}_discovery.txt")

    # ---- clonal redundancy --------------------------------------
    clusters = dereplicate(tree, a.derep_thresh)
    memb = {re.sub(r"\.(fna|fa|fasta)(\.gz)?$", "", g): i
            for i, c in enumerate(clusters) for g in c}
    meta["clone_cluster"] = meta.genome.map(memb)
    meta["clone_size"] = meta.clone_cluster.map(meta.groupby("clone_cluster").size())
    meta["weight"] = 1.0 / meta.clone_size

    n_multi = (meta.groupby("clone_cluster").size() > 1).sum()
    print(f"\nclonal clusters (threshold {a.derep_thresh}): "
          f"{meta.clone_cluster.nunique()} total, {n_multi} with >1 member")
    by_group = meta.groupby("fastbaps_group")["clone_size"].mean()
    print("mean clonal cluster size per group (focal first):")
    print(by_group.reindex([focal_label] +
                           [g for g in by_group.index if g != focal_label]).head(15))
    if by_group.max() / max(by_group.min(), 1e-9) > 2:
        print("\nWARNING: clonal redundancy is uneven across groups. An apparent "
              "IS excess in the focal group can be produced entirely by this. "
              "Use the weights in script 02, and repeat on one representative "
              "per cluster as a sensitivity check.")

    meta[["genome", "fastbaps_group", "clone_cluster", "weight"]].to_csv(
        f"{a.out_prefix}_weights.tsv", sep="\t", index=False)

    # ---- scoring panel ------------------------------------------
    focal = meta[meta.fastbaps_group == focal_label]
    nonfocal = meta[meta.fastbaps_group != focal_label]

    keep = [focal]   # every focal genome: the SFS is computed here
    if len(nonfocal) > a.max_nonfocal:
        reps = nonfocal.sort_values("weight", ascending=False) \
                       .drop_duplicates("clone_cluster")
        per_group = max(5, a.max_nonfocal // max(reps.fastbaps_group.nunique(), 1))
        samp = (reps.groupby("fastbaps_group", group_keys=False)
                    .apply(lambda x: x.sample(min(len(x), per_group),
                                              random_state=1)))
        keep.append(samp)
        print(f"\nscoring panel: all {len(focal)} focal + {len(samp)} non-focal "
              f"(stratified, up to {per_group} per group)")
    else:
        keep.append(nonfocal)
        print(f"\nscoring panel: all {len(focal)} focal + {len(nonfocal)} non-focal")

    panel = pd.concat(keep)
    panel.genome.to_csv(f"{a.out_prefix}_scoring.txt", index=False, header=False)

    thin = (panel[panel.fastbaps_group != focal_label]
            .groupby("fastbaps_group").size())
    thin = thin[thin < 5]
    if len(thin):
        print("\nWARNING: these groups have <5 genomes in the scoring panel; "
              "insertions will often be unpolarisable against them:")
        print(thin)

    # ---- what script 02 will additionally need ------------------
    need = {"n_contigs", "N50", "total_len"}
    have = {c.lower() for c in pd.read_csv(a.meta, sep=sniff_sep(a.meta),
                                           nrows=0).columns}
    missing = {c for c in need if c.lower() not in have}
    if missing:
        print(f"\nnote: {sorted(missing)} are not in this metadata file. "
              "Script 02 needs them as covariates so that assembly "
              "fragmentation cannot masquerade as an IS-load difference. "
              "Generate them with:\n"
              "  bash 00a_prep_inputs.sh genomes <your_fastbaps_file>")


if __name__ == "__main__":
    main()
