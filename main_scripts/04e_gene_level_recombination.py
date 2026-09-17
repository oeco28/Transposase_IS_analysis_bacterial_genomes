#!/usr/bin/env python3
"""
04e_gene_level_recombination.py

Replaces test C3 in 04c when gubbins was run on a CONCATENATED CORE GENE
alignment rather than a whole-genome reference.

The problem with coordinates: a concatenated core alignment contains only
core genes, no intergenic sequence. IS elements almost never sit inside
core genes, so IS positions cannot be placed on that axis at all.

The fix: join on GENES instead of base positions.
  1. Alignment position -> core gene, via the partition file your
     concatenation tool wrote (RAxML/IQ-TREE partitions, or any table of
     gene, start, end).
  2. Recombination block -> the set of core genes it covers.
  3. IS copy -> its neighbouring core genes, from the per-genome GFFs and
     the pangenome clustering (real genome coordinates, where IS elements
     actually live).
  4. Test: are IS-adjacent core genes over-represented among genes hit by
     recombination, relative to core genes that are not IS-adjacent?

Mechanistic reading: IS copies are dispersed homologous repeats. If they
template ectopic exchange, recombination tracts should terminate near
them, so their FLANKING core genes should be recombination-rich. This is
a stronger prediction than blocks simply containing IS elements.

Inputs
  --partitions   gene -> alignment coordinates (partition file or TSV)
  --recomb-gff   one or more gubbins recombination_predictions.gff
  --gff          directory of per-genome GFFs (real coordinates)
  --calls        is_out/IS_calls.tsv
  --pa           gene_presence_absence.csv (to map locus tags to clusters)
  --meta         metadata.tsv

Usage:
  python 04e_gene_level_recombination.py --partitions core_genes.partitions \
      --gubbins-dir gubbins --gff gff --calls is_out/IS_calls.tsv \
      --pa panaroo_out/gene_presence_absence.csv --meta metadata.tsv --focal 12
"""

import argparse
import glob
import os
import re
import sys

import numpy as np
import pandas as pd
from scipy import stats


def parse_partitions(path):
    """
    Accepts RAxML/IQ-TREE partition lines:
        DNA, gene_0001 = 1-1071
        charset gene_0001 = 1-1071;
    or a TSV/CSV with columns like gene,start,end.
    """
    rows = []
    with open(path) as fh:
        txt = fh.read()
    pat = re.compile(r"([\w.\-]+)\s*=\s*(\d+)\s*-\s*(\d+)")
    for m in pat.finditer(txt):
        rows.append((m.group(1), int(m.group(2)), int(m.group(3))))
    if rows:
        return pd.DataFrame(rows, columns=["gene", "start", "end"])

    sep = "\t" if path.endswith((".tsv", ".txt")) else ","
    t = pd.read_csv(path, sep=sep)
    cols = {c.lower(): c for c in t.columns}
    g = next((cols[c] for c in cols if "gene" in c or "clust" in c), None)
    s = next((cols[c] for c in cols if c.startswith("start") or c == "from"), None)
    e = next((cols[c] for c in cols if c.startswith("end") or c == "to"), None)
    if not all([g, s, e]):
        sys.exit(f"could not parse {path}. Provide gene/start/end columns, or a "
                 f"RAxML-style partition file. Columns found: {list(t.columns)}")
    return t[[g, s, e]].rename(columns={g: "gene", s: "start", e: "end"})


def parse_recomb_gff(path):
    rows = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 5:
                continue
            try:
                rows.append((int(f[3]), int(f[4])))
            except ValueError:
                continue
    return pd.DataFrame(rows, columns=["start", "end"])


def genes_in_blocks(blocks, parts):
    """Count how many recombination blocks overlap each core gene."""
    starts = parts.start.values
    ends = parts.end.values
    hits = np.zeros(len(parts), dtype=int)
    for s, e in zip(blocks.start.values, blocks.end.values):
        hits += (starts <= e) & (ends >= s)
    out = parts.copy()
    out["n_blocks"] = hits
    return out


def parse_gff_cds(path):
    rows = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("##FASTA"):
                break
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "CDS":
                continue
            m = re.search(r"ID=([^;]+)", f[8])
            rows.append((f[0], int(f[3]), int(f[4]), m.group(1) if m else None))
    return pd.DataFrame(rows, columns=["contig", "start", "end", "locus_tag"])


def is_adjacent_clusters(calls, gff_dir, locus2cluster, genomes, window):
    """Core-gene clusters within `window` bp of an IS copy, in any genome."""
    near = set()
    seen = 0
    for g in genomes:
        p = os.path.join(gff_dir, f"{g}.gff")
        if not os.path.exists(p):
            continue
        cds = parse_gff_cds(p)
        if cds.empty:
            continue
        ist = calls[calls.genome == g]
        if ist.empty:
            continue
        seen += 1
        for contig, sub in cds.groupby("contig"):
            ic = ist[ist.contig == contig]
            if ic.empty:
                continue
            s_arr, e_arr = ic.start.values, ic.end.values
            for _, r in sub.iterrows():
                if ((s_arr < r.end + window) & (e_arr > r.start - window)).any():
                    cl = locus2cluster.get(r.locus_tag)
                    if cl:
                        near.add(cl)
    return near, seen


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--partitions", required=True)
    ap.add_argument("--gubbins-dir", default=None)
    ap.add_argument("--recomb-gff", nargs="*", default=None)
    ap.add_argument("--gff", default="gff")
    ap.add_argument("--calls", default="is_out/IS_calls.tsv")
    ap.add_argument("--pa", default="gene_presence_absence.csv")
    ap.add_argument("--meta", default="metadata.tsv")
    ap.add_argument("--focal", required=True)
    ap.add_argument("--window", type=int, default=10000,
                    help="bp around an IS copy that counts as adjacent")
    ap.add_argument("--out", default="is_out")
    a = ap.parse_args()

    os.makedirs(a.out, exist_ok=True)
    parts = parse_partitions(a.partitions)
    print(f"{len(parts)} core genes in the concatenated alignment "
          f"(length {int(parts.end.max())} bp)")

    # ---- recombination per core gene, per group -----------------
    gffs = {}
    if a.recomb_gff:
        for p in a.recomb_gff:
            gffs[os.path.basename(os.path.dirname(p)) or os.path.basename(p)] = p
    elif a.gubbins_dir:
        for d in sorted(glob.glob(os.path.join(a.gubbins_dir, "*"))):
            hit = glob.glob(os.path.join(d, "*recombination_predictions.gff"))
            if hit:
                m = re.search(r"(\d+)\s*$", os.path.basename(d))
                gffs[m.group(1) if m else os.path.basename(d)] = hit[0]
    if not gffs:
        sys.exit("no recombination GFFs found; give --gubbins-dir or --recomb-gff")
    print(f"recombination predictions for {len(gffs)} groups: {sorted(gffs)}")

    per_gene = {}
    for grp, p in gffs.items():
        b = parse_recomb_gff(p)
        if b.empty:
            continue
        gg = genes_in_blocks(b, parts)
        per_gene[grp] = gg.set_index("gene")["n_blocks"]
        print(f"  group {grp}: {len(b)} blocks, "
              f"{(gg.n_blocks > 0).sum()} of {len(gg)} core genes hit")
    G = pd.DataFrame(per_gene)
    G.to_csv(os.path.join(a.out, "recombination_per_core_gene.tsv"), sep="\t")

    # ---- which core genes are IS-adjacent -----------------------
    meta = pd.read_csv(a.meta, sep="\t", dtype={"fastbaps_group": str},
                       low_memory=False)
    meta["genome"] = meta.genome.astype(str)
    meta["fastbaps_group"] = meta.fastbaps_group.astype(str).str.strip()
    focal = str(a.focal).strip()
    focal_genomes = meta.loc[meta.fastbaps_group == focal, "genome"].tolist()

    pa = pd.read_csv(a.pa, low_memory=False, dtype=str)
    cluster_col = pa.columns[0]
    genome_cols = [c for c in pa.columns
                   if re.sub(r"\.(gff|fna|fa|fasta)$", "", str(c)) in set(meta.genome)
                   or str(c) in set(meta.genome)]
    print(f"pangenome: {len(pa)} clusters, {len(genome_cols)} genome columns")
    locus2cluster = {}
    for c in genome_cols:
        for cl, v in zip(pa[cluster_col], pa[c]):
            if isinstance(v, str):
                for lt in re.split(r"[;\t ]+", v):
                    if lt:
                        locus2cluster[lt] = cl

    calls = pd.read_csv(a.calls, sep="\t", low_memory=False)
    calls["genome"] = calls.genome.astype(str)

    print(f"focal group '{focal}': {len(focal_genomes)} genomes")
    if not focal_genomes:
        sys.exit(f"no genomes in group '{focal}'; groups present: "
                 f"{sorted(set(meta.fastbaps_group))}")
    if not locus2cluster:
        sys.exit("no locus tags could be read from the presence/absence table; "
                 "run 04a_check_pangenome_inputs.py to diagnose")

    near, n_used = is_adjacent_clusters(calls, a.gff, locus2cluster,
                                        focal_genomes, a.window)
    print(f"scanned {n_used} focal genomes; {len(near)} clusters lie within "
          f"{a.window} bp of an IS copy")
    if n_used == 0:
        sys.exit(f"none of the {len(focal_genomes)} focal genomes had both a GFF "
                 f"in {a.gff}/ and IS calls.\n"
                 f"  GFFs present: {len(glob.glob(os.path.join(a.gff, '*.gff')))}\n"
                 f"  genomes in IS calls: {calls.genome.nunique()}\n"
                 f"  Run 04a_check_pangenome_inputs.py to find the mismatch.")

    core_genes = set(parts.gene)
    near_core = near & core_genes
    print(f"of those, {len(near_core)} are core genes in the alignment")
    if len(near_core) < 10:
        print("\nToo few IS-adjacent core genes for a test. This is the expected "
              "result if\nIS elements sit almost entirely in accessory regions, "
              "which is itself worth\nreporting: it means transposition and core "
              "recombination are acting on\nlargely disjoint parts of the genome.")
        return

    # ---- the test -----------------------------------------------
    print("\n=== recombination in IS-adjacent vs other core genes ===")
    rows = []
    for grp in G.columns:
        v = G[grp].dropna()
        adj = v[v.index.isin(near_core)]
        oth = v[~v.index.isin(near_core)]
        if len(adj) < 5 or len(oth) < 5:
            continue
        u = stats.mannwhitneyu(adj, oth, alternative="two-sided")
        rows.append({"group": grp, "n_adjacent": len(adj), "n_other": len(oth),
                     "median_adj": float(adj.median()),
                     "median_other": float(oth.median()),
                     "mean_adj": float(adj.mean()),
                     "mean_other": float(oth.mean()),
                     "p": u.pvalue})
    if not rows:
        print("  not enough genes in either category")
        return
    res = pd.DataFrame(rows)
    res["enriched"] = res.mean_adj > res.mean_other
    res = res.sort_values("p")
    print(res.round(4).to_string(index=False))
    res.to_csv(os.path.join(a.out, "IS_adjacent_core_gene_recombination.tsv"),
               sep="\t", index=False)

    n_up = int(res.enriched.sum())
    sign_p = stats.binomtest(n_up, len(res), 0.5).pvalue
    print(f"\n  IS-adjacent core genes recombine more in {n_up}/{len(res)} "
          f"groups; sign test p = {sign_p:.4g}")
    print("  Consistency of direction across independent lineages is the result "
          "to report.\n  A focal-only effect says the association is specific to "
          f"group {focal};\n  a species-wide effect says IS elements mark "
          "recombination margins generally.")

    print("\nnote: core-genome recombination only, and IS adjacency was scored "
          "in focal\ngenomes. Both belong in the methods.")


if __name__ == "__main__":
    main()
