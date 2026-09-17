#!/usr/bin/env python3
"""
00b_merge_ncbi_metadata.py

Turns an NCBI Datasets assembly table plus your fastbaps output into the
canonical metadata.tsv the rest of the pipeline expects:

    genome  fastbaps_group  n_contigs  N50  total_len  assembly_level
    checkm_completeness  checkm_contamination

Why a dedicated script rather than hand-recoding the columns:

  * The NCBI table carries assembly-quality fields that are not
    decoration here. They are the covariates that decide whether your
    IS-load result is real. Recoding by hand tends to drop them.

  * 'Assembly Level' is the single most dangerous confounder in this
    project. Complete genomes recover far more IS copies than draft
    assemblies, because short-read contigs break at repeats and collapse
    identical copies. If assembly level is unevenly distributed across
    fastbaps groups (very common: one lineage gets long-read sequenced
    by one consortium), that alone can manufacture your entire result.
    This script quantifies that imbalance so you find out now.

  * CheckM contamination directly inflates IS counts, because a
    contaminating contig brings its own IS elements with it. Anything
    above ~5% contamination should be dropped, not modelled.

  * NCBI's table has no total sequence length in the default column set,
    so genome size is computed from the assemblies with seqkit.

Usage:
  python 00b_merge_ncbi_metadata.py \
      --ncbi ncbi_assembly_table.tsv \
      --fastbaps fastbaps_clusters.csv \
      --asm assemblies \
      --out metadata.tsv
"""

import argparse
import io
import os
import re
import subprocess
import sys

import pandas as pd

# NCBI Datasets column name -> canonical name
NCBI_MAP = {
    "Assembly Accession": "genome",
    "Assembly Name": "assembly_name",
    "Organism Infraspecific Names Strain": "strain",
    "Organism Infraspecific Names Isolate": "isolate",
    "Assembly Level": "assembly_level",
    "Assembly Stats Scaffold N50": "N50",
    "Assembly Stats Number of Scaffolds": "n_contigs",
    "Assembly Sequencing Tech": "seq_tech",
    "CheckM completeness": "checkm_completeness",
    "CheckM contamination": "checkm_contamination",
    "Assembly Release Date": "release_date",
    "WGS project accession": "wgs_accession",
}


def sniff(path):
    if path.endswith((".tsv", ".txt", ".tab")):
        return "\t"
    if path.endswith(".csv"):
        return ","
    with open(path) as fh:
        h = fh.readline()
    return "\t" if h.count("\t") > h.count(",") else ","


def norm_id(s):
    """Strip file extensions and the NCBI _ASM..._genomic suffix."""
    return (s.astype(str).str.strip()
             .str.replace(r"_genomic$", "", regex=True)
             .str.replace(r"\.(fna|fa|fasta|gbff)(\.gz)?$", "", regex=True)
             .str.replace(r"^(GC[AF]_\d+\.\d+)_.*$", r"\1", regex=True))


def seqkit_lengths(asm_dir):
    """Total assembly length per genome; NCBI's default columns omit it."""
    files = [f for f in os.listdir(asm_dir) if f.endswith(".fna")]
    if not files:
        return None
    try:
        out = subprocess.run(
            ["seqkit", "stats", "-a", "-T"] + [os.path.join(asm_dir, f) for f in files],
            capture_output=True, text=True, check=True).stdout
    except FileNotFoundError:
        print("  (seqkit is not on PATH: conda activate ismobile)")
        return None
    except subprocess.CalledProcessError as e:
        print(f"  (seqkit failed: {e.stderr.strip()[:200]})")
        return None
    st = pd.read_csv(io.StringIO(out), sep="\t")
    st["genome"] = st["file"].apply(lambda p: os.path.basename(p).rsplit(".", 1)[0])
    return st[["genome", "sum_len", "num_seqs", "N50"]].rename(
        columns={"sum_len": "total_len_obs", "num_seqs": "n_contigs_obs",
                 "N50": "N50_obs"})


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ncbi", required=True)
    ap.add_argument("--fastbaps", required=True)
    ap.add_argument("--asm", default="assemblies",
                    help="directory of <genome>.fna, for total_len")
    ap.add_argument("--out", default="metadata.tsv")
    ap.add_argument("--fb-id-col", default=None)
    ap.add_argument("--fb-group-col", default=None)
    ap.add_argument("--max-contamination", type=float, default=5.0)
    ap.add_argument("--min-completeness", type=float, default=90.0)
    ap.add_argument("--no-filter", action="store_true",
                    help="report QC failures but keep them in the output")
    a = ap.parse_args()

    # ---- NCBI table --------------------------------------------
    nc = pd.read_csv(a.ncbi, sep=sniff(a.ncbi), dtype=str, low_memory=False)
    present = {k: v for k, v in NCBI_MAP.items() if k in nc.columns}
    missing = [k for k in ("Assembly Accession", "Assembly Level")
               if k not in nc.columns]
    if missing:
        sys.exit(f"NCBI table is missing required columns: {missing}\n"
                 f"Columns present: {list(nc.columns)}")
    nc = nc.rename(columns=present)[list(present.values())]
    nc["genome"] = norm_id(nc["genome"])
    for c in ("N50", "n_contigs", "checkm_completeness", "checkm_contamination"):
        if c in nc.columns:
            nc[c] = pd.to_numeric(nc[c], errors="coerce")
    print(f"NCBI table: {len(nc)} assemblies, "
          f"{len(present)} of {len(NCBI_MAP)} known columns mapped")

    # ---- fastbaps ----------------------------------------------
    fb = pd.read_csv(a.fastbaps, sep=sniff(a.fastbaps), dtype=str)

    # R's write.csv emits an unnamed rowname column, which pandas calls
    # "Unnamed: 0". It contains "name", so a naive regex picks it as the
    # genome ID and the join silently matches nothing. Drop those first.
    junk = [c for c in fb.columns
            if re.match(r"^(unnamed|x)(:?\s*\d*)?$", str(c).strip(), re.I)
            or str(c).strip() == ""]
    if junk:
        print(f"fastbaps: ignoring index-like column(s) {junk}")
        fb = fb.drop(columns=junk)

    id_cands = [c for c in fb.columns
                if re.search(r"isolate|genome|sample|taxon|accession|strain|^id$|name",
                             str(c), re.I)]
    grp_cands = [c for c in fb.columns
                 if re.search(r"fastbaps|baps|cluster|group|partition|lineage",
                              str(c), re.I)]

    # Validate anything passed explicitly, before pandas turns it into a
    # bare KeyError. Note that --focal belongs to 00_build_panel.py, not
    # here; this script does not need to know the focal group.
    for flag, val in (("--fb-id-col", a.fb_id_col), ("--fb-group-col", a.fb_group_col)):
        if val is not None and val not in fb.columns:
            hint = ""
            if " " in val and val.split()[0] in fb.columns:
                hint = (f"\n\nIt looks like '{val.split()[0]}' is the column and "
                        f"'{' '.join(val.split()[1:])}' was meant as the focal group. "
                        f"This script has no --focal option; the focal group is given "
                        f"to 00_build_panel.py in the next step.")
            sys.exit(f"\n{flag} '{val}' is not a column in {a.fastbaps}.\n"
                     f"Columns present: {list(fb.columns)}{hint}")

    idc = a.fb_id_col or (id_cands[0] if id_cands else None)
    grc = a.fb_group_col
    if grc is None:
        if len(grp_cands) > 1:
            sys.exit("\nfastbaps has several clustering columns:\n" +
                     "\n".join(f"    {c}" for c in grp_cands) +
                     "\n\nThese are nested BAPS levels, not alternatives: level 1 is a "
                     "coarse split and each higher level subdivides it. Picking the "
                     "wrong one silently changes what 'the focal group' means.\n"
                     "Choose the level your phylogenetic comparison was done at:\n"
                     f"    --fb-group-col '{grp_cands[-1]}'")
        grc = grp_cands[0] if grp_cands else None

    if idc is None or grc is None:
        sys.exit(f"Could not detect ID/group columns in {a.fastbaps}.\n"
                 f"Columns: {list(fb.columns)}\n"
                 f"Use --fb-id-col and --fb-group-col.")
    if len(id_cands) > 1 and a.fb_id_col is None:
        print(f"note: several candidate ID columns {id_cands}; using '{idc}'")
    print(f"fastbaps: ID column '{idc}', group column '{grc}'")
    fb = fb[[idc, grc]].rename(columns={idc: "genome", grc: "fastbaps_group"})
    fb["genome"] = norm_id(fb["genome"])
    fb["fastbaps_group"] = fb["fastbaps_group"].astype(str).str.strip()

    # ---- join ---------------------------------------------------
    m = nc.merge(fb, on="genome", how="outer", indicator=True)
    n_ncbi_only = (m._merge == "left_only").sum()
    n_fb_only = (m._merge == "right_only").sum()
    if n_ncbi_only:
        print(f"WARNING: {n_ncbi_only} NCBI assemblies have no fastbaps group, "
              f"e.g. {m.loc[m._merge=='left_only','genome'].head(3).tolist()}")
    if n_fb_only:
        print(f"WARNING: {n_fb_only} fastbaps entries have no NCBI row, "
              f"e.g. {m.loc[m._merge=='right_only','genome'].head(3).tolist()}")
    if n_ncbi_only or n_fb_only:
        print("  This is usually an accession-format mismatch (GCA vs GCF, or "
              "version suffix), not missing data. Check before continuing.")
    m = m[m._merge == "both"].drop(columns="_merge")
    print(f"joined: {len(m)} genomes")

    # ---- genome length from the assemblies ----------------------
    lens = seqkit_lengths(a.asm) if os.path.isdir(a.asm) else None
    if lens is not None:
        m = m.merge(lens, on="genome", how="left")
        m["total_len"] = m["total_len_obs"]
        # prefer observed contig counts: NCBI counts scaffolds, and the
        # analysis cares about the contigs ISEScan actually sees
        m["n_contigs"] = m["n_contigs_obs"].fillna(m.get("n_contigs"))
        m["N50"] = m["N50_obs"].fillna(m.get("N50"))
        m = m.drop(columns=[c for c in ("total_len_obs", "n_contigs_obs", "N50_obs")
                            if c in m.columns])
        print(f"computed total_len from {a.asm}/ for {m.total_len.notna().sum()} genomes")
    else:
        m["total_len"] = pd.NA
        print(f"WARNING: could not read {a.asm}/ with seqkit, so total_len is empty.\n"
              "  Script 02 uses log(total_len) as the model offset and will fail.\n"
              "  Run 00a_prep_inputs.sh first, or activate the ismobile env.")

    # ---- QC: contamination --------------------------------------
    if "checkm_contamination" in m.columns:
        bad = m[(m.checkm_contamination > a.max_contamination) |
                (m.checkm_completeness < a.min_completeness)]
        print(f"\nQC: {len(bad)} genomes exceed {a.max_contamination}% contamination "
              f"or fall below {a.min_completeness}% completeness")
        if len(bad):
            print("  Contaminating contigs carry their own IS elements, so these "
                  "inflate IS counts directly. Drop them rather than modelling them.")
            if len(bad) and "fastbaps_group" in bad.columns:
                print("  by group:")
                print(bad.groupby("fastbaps_group").size().to_string())
        if not a.no_filter:
            m = m.drop(bad.index)
            print(f"  removed; {len(m)} genomes retained (use --no-filter to keep)")

    # ---- the confounder that matters ----------------------------
    print("\n=== assembly level composition per fastbaps group ===")
    lvl = (m.groupby(["fastbaps_group", "assembly_level"]).size()
             .unstack(fill_value=0))
    frac_complete = (lvl.get("Complete Genome", 0) / lvl.sum(axis=1)).sort_values()
    print(lvl.to_string())
    print("\nfraction of complete genomes per group:")
    print(frac_complete.round(3).to_string())
    spread = frac_complete.max() - frac_complete.min()
    if spread > 0.15:
        print(f"\n*** WARNING: complete-genome fraction ranges over {spread:.2f} "
              "across groups.\n"
              "    Complete assemblies recover substantially more IS copies than\n"
              "    draft ones, because short-read contigs break at repeats and\n"
              "    collapse identical copies. If your focal group happens to be\n"
              "    better assembled, that alone can produce an apparent IS excess.\n"
              "    Mitigations, in order of strength:\n"
              "      1. restrict the primary analysis to one assembly level\n"
              "      2. add assembly_level as a covariate in script 02\n"
              "      3. report the result stratified by assembly level")

    if "seq_tech" in m.columns:
        tech = m.groupby("fastbaps_group")["seq_tech"].apply(
            lambda x: x.astype(str).str.contains("nanopore|pacbio|smrt", case=False).mean())
        if tech.max() - tech.min() > 0.15:
            print("\nWARNING: long-read sequencing is also unevenly distributed "
                  "across groups:")
            print(tech.round(3).to_string())

    # ---- write --------------------------------------------------
    cols = ["genome", "fastbaps_group", "n_contigs", "N50", "total_len",
            "assembly_level", "checkm_completeness", "checkm_contamination",
            "seq_tech", "strain", "isolate"]
    out = m[[c for c in cols if c in m.columns]]
    out.to_csv(a.out, sep="\t", index=False)
    print(f"\nwrote {a.out}: {len(out)} genomes, {out.fastbaps_group.nunique()} groups")
    print(out.fastbaps_group.value_counts().head(20).to_string())


if __name__ == "__main__":
    main()
