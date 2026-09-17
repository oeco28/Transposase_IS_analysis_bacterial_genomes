#!/usr/bin/env python3
"""
03_insertion_site_polymorphism.py

Turns IS calls into a *population-genetic* dataset: a genome x
insertion-site occupancy matrix, scored as occupied / empty / unknown.

Why this is the analysis that makes the story:
  - Presence/absence of an accessory GENE is ambiguous (gain by
    transposition, gain by homologous recombination, or loss).
  - Presence/absence of an IS at a defined nucleotide SITE, with the
    empty allele observed in sister genomes, is direct evidence that a
    transposition event happened on a specific branch.
  - Once you have the matrix you can do everything you would normally
    do with SNPs: frequency spectrum, homoplasy, ancestral
    reconstruction, Fst, GWAS.

Method (assemblies only, no reads needed):
  For each IS copy, cut L bp of left and right flank. Build two probes:
     occupied probe = left flank + first 150 bp of the IS
     empty probe    = left flank + right flank (junction, no IS)
  BLAST both against every genome.
     full-length hit to empty probe  -> EMPTY   (ancestral, no insertion)
     full-length hit to occupied probe -> OCCUPIED
     neither -> UNKNOWN (region absent, or contig break: excluded)
  Cluster probes whose flanks are homologous into one site.

If you also have short reads, run panISa or ISMapper in parallel; they
detect insertions present in only part of the population within a
sample, which is an even stronger argument for ongoing activity.
"""

import argparse, os, subprocess, tempfile, sys
from collections import defaultdict
import pandas as pd
import numpy as np
from Bio import SeqIO

FLANK = 400
IS_TAG = 150
MIN_COV = 0.92
MIN_ID = 95.0


def load_calls(path):
    d = pd.read_csv(path, sep="\t")
    d = d[d["complete"]] if "complete" in d.columns else d
    return d


def build_probes(calls, asm_dir, out_fa):
    """One probe pair per IS copy, skipping copies too close to a contig end."""
    recs = []
    for g, sub in calls.groupby("genome"):
        idx = SeqIO.index(os.path.join(asm_dir, f"{g}.fna"), "fasta")
        for _, r in sub.iterrows():
            seq = idx[r["contig"]].seq
            s, e = int(r["start"]) - 1, int(r["end"])
            if s < FLANK or e + FLANK > len(seq):
                continue                      # contig edge: not scorable
            L, R = seq[s - FLANK:s], seq[e:e + FLANK]
            pid = f"{g}::{r['contig']}::{s}::{r['family']}"
            recs.append((f"{pid}::EMPTY", str(L) + str(R)))
            recs.append((f"{pid}::OCCUPIED", str(L) + str(seq[s:s + IS_TAG])))
    with open(out_fa, "w") as fh:
        for n, s in recs:
            fh.write(f">{n}\n{s}\n")
    return len(recs) // 2


def blast_all(probe_fa, asm_dir, genomes, threads, workdir):
    """BLAST the whole probe set against each genome; keep full-length hits."""
    hits = []
    for g in genomes:
        db = os.path.join(workdir, g)
        subprocess.run(["makeblastdb", "-in", os.path.join(asm_dir, f"{g}.fna"),
                        "-dbtype", "nucl", "-out", db],
                       check=True, stdout=subprocess.DEVNULL)
        res = subprocess.run(
            ["blastn", "-query", probe_fa, "-db", db, "-num_threads", str(threads),
             "-perc_identity", str(MIN_ID), "-max_target_seqs", "5",
             "-outfmt", "6 qseqid sseqid pident length qlen sstart send"],
            check=True, capture_output=True, text=True).stdout
        for line in res.strip().split("\n"):
            if not line:
                continue
            q, s, pid, alen, qlen, ss, se = line.split("\t")
            if int(alen) >= MIN_COV * int(qlen):
                hits.append((g, q, float(pid)))
    return pd.DataFrame(hits, columns=["genome", "probe", "pident"])


def score_matrix(hits, genomes):
    """occupied=1, empty=0, unknown=NaN."""
    hits[["site", "allele"]] = hits["probe"].str.rsplit("::", n=1, expand=True)
    best = (hits.sort_values("pident", ascending=False)
                .drop_duplicates(["genome", "site", "allele"]))
    piv = best.pivot_table(index="site", columns=["genome", "allele"],
                           values="pident", aggfunc="max")
    sites = piv.index
    mat = pd.DataFrame(np.nan, index=sites, columns=genomes, dtype=float)
    for g in genomes:
        occ = piv[(g, "OCCUPIED")] if (g, "OCCUPIED") in piv.columns else pd.Series(np.nan, index=sites)
        emp = piv[(g, "EMPTY")] if (g, "EMPTY") in piv.columns else pd.Series(np.nan, index=sites)
        mat.loc[occ.notna() & emp.isna(), g] = 1.0
        mat.loc[emp.notna() & occ.isna(), g] = 0.0
        # both hit: tandem/duplicated site, leave unknown rather than guess
    return mat


def collapse_sites(mat, jaccard=0.95):
    """Different genomes rediscover the same insertion; merge identical rows."""
    key = mat.fillna(-1).astype(int).astype(str).agg("".join, axis=1)
    keep = mat.groupby(key.values).first()
    keep.index = [f"site{i:05d}" for i in range(len(keep))]
    return keep


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--calls", default="is_out/IS_calls.tsv")
    ap.add_argument("--asm", default="assemblies")
    ap.add_argument("--meta", default="metadata.tsv")
    ap.add_argument("--out", default="is_out/insertion_matrix.tsv")
    ap.add_argument("--threads", type=int, default=8)
    a = ap.parse_args()

    calls = load_calls(a.calls)
    meta = pd.read_csv(a.meta, sep="\t")
    genomes = sorted(meta["genome"].astype(str))

    with tempfile.TemporaryDirectory() as wd:
        probe_fa = os.path.join(wd, "probes.fna")
        n = build_probes(calls, a.asm, probe_fa)
        print(f"{n} scorable IS copies -> {2*n} probes", file=sys.stderr)
        hits = blast_all(probe_fa, a.asm, genomes, a.threads, wd)

    mat = collapse_sites(score_matrix(hits, genomes))
    mat.to_csv(a.out, sep="\t")

    # ---- frequency spectrum -------------------------------------
    # Neutral transposition-selection balance predicts an SFS strongly
    # skewed to singletons (most insertions are deleterious or neutral
    # and young). Insertions segregating at high frequency, or fixed in
    # a subclade, are your adaptive candidates. Compare the IS SFS to
    # the SFS of synonymous core SNPs from the same genomes: a
    # significant excess of high-frequency insertions relative to that
    # neutral reference is the population-genetic claim you want to make.
    focal = set(meta.loc[meta.fastbaps_group == "G4", "genome"].astype(str))
    sub = mat[[c for c in mat.columns if c in focal]]
    freq = sub.mean(axis=1, skipna=True)
    n_scored = sub.notna().sum(axis=1)
    out = pd.DataFrame({"freq": freq, "n_scored": n_scored})
    out = out[out.n_scored >= 0.5 * len(focal)]
    out.to_csv("is_out/insertion_frequencies.tsv", sep="\t")
    print(out["freq"].describe())


if __name__ == "__main__":
    main()
