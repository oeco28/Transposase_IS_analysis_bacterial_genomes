#!/usr/bin/env bash
# ==================================================================
# 00a_prep_inputs.sh
#
# Two jobs, both of which the later scripts assume are already done:
#
#   1. Normalise assembly filenames. Script 01 globs *.fna and derives
#      the genome ID from the basename. Your files are probably .fasta
#      or .fa, and their basenames may not match the IDs in your
#      fastbaps output or your tree tips. This creates a directory of
#      SYMLINKS with consistent names, leaving your genomes/ untouched.
#
#   2. Build metadata.tsv, which every downstream script reads:
#         genome  fastbaps_group  n_contigs  N50  total_len
#      The assembly stats are not optional decoration. They are the
#      covariates that stop assembly fragmentation masquerading as an
#      IS-load difference between groups.
#
# Usage:
#   bash 00a_prep_inputs.sh genomes fastbaps_clusters.csv
#   bash 00a_prep_inputs.sh genomes fastbaps_clusters.csv assemblies metadata.tsv
#
# The fastbaps file is expected to have a genome-ID column and a cluster
# column; the script tries to detect them and tells you if it cannot.
# ==================================================================
set -euo pipefail

GENOMES=${1:-genomes}
FASTBAPS=${2:-fastbaps_clusters.csv}
LINKDIR=${3:-assemblies}
META=${4:-metadata.tsv}

command -v seqkit >/dev/null || { echo "seqkit not found: conda activate ismobile"; exit 1; }
[[ -d "$GENOMES" ]] || { echo "no such directory: $GENOMES"; exit 1; }
[[ -f "$FASTBAPS" ]] || { echo "no such file: $FASTBAPS"; exit 1; }

# ---- 1. symlink farm with .fna extensions ------------------------
mkdir -p "$LINKDIR"
n=0
shopt -s nullglob
for f in "$GENOMES"/*.{fna,fa,fasta,fna.gz,fa.gz,fasta.gz}; do
  base=$(basename "$f")
  # strip one or two extensions to get the genome ID
  id="${base%.gz}"; id="${id%.*}"
  if [[ "$f" == *.gz ]]; then
    # ISEScan will not read gzipped input; decompress these for real
    [[ -f "$LINKDIR/$id.fna" ]] || zcat "$f" > "$LINKDIR/$id.fna"
  else
    ln -sf "$(readlink -f "$f")" "$LINKDIR/$id.fna"
  fi
  n=$((n+1))
done
shopt -u nullglob
echo "linked/prepared $n assemblies into $LINKDIR/"
[[ $n -gt 0 ]] || { echo "found no assemblies in $GENOMES"; exit 1; }

# ---- 2. assembly statistics --------------------------------------
echo "computing assembly stats (this reads every genome once)"
seqkit stats -a -T "$LINKDIR"/*.fna > .seqkit_stats.tsv

# ---- 3. join to fastbaps groups ----------------------------------
python3 - "$FASTBAPS" "$META" <<'PY'
import sys, os, re
import pandas as pd

fastbaps_path, out = sys.argv[1], sys.argv[2]

st = pd.read_csv(".seqkit_stats.tsv", sep="\t")
st["genome"] = st["file"].apply(lambda p: os.path.basename(p).rsplit(".", 1)[0])
st = st.rename(columns={"num_seqs": "n_contigs", "sum_len": "total_len", "N50": "N50"})
st = st[["genome", "n_contigs", "total_len", "N50"]]

sep = "\t" if fastbaps_path.endswith((".tsv", ".txt")) else ","
fb = pd.read_csv(fastbaps_path, sep=sep)

# fastbaps output column names vary by how it was run; detect rather than assume
id_col = next((c for c in fb.columns
               if re.search(r"isolate|genome|sample|taxon|id|name", c, re.I)), None)
grp_col = next((c for c in fb.columns
                if re.search(r"clust|group|baps|level|partition", c, re.I)), None)
if id_col is None or grp_col is None:
    sys.exit(f"could not detect ID/cluster columns in {fastbaps_path}. "
             f"Columns present: {list(fb.columns)}. "
             f"Rename them to 'genome' and 'fastbaps_group' and rerun.")
print(f"using '{id_col}' as genome ID and '{grp_col}' as fastbaps group")

fb = fb[[id_col, grp_col]].rename(columns={id_col: "genome", grp_col: "fastbaps_group"})
fb["genome"] = fb["genome"].astype(str).str.replace(r"\.(fna|fa|fasta)$", "", regex=True)
fb["fastbaps_group"] = fb["fastbaps_group"].astype(str)

m = st.merge(fb, on="genome", how="outer", indicator=True)

only_asm = m[m._merge == "left_only"]["genome"].tolist()
only_fb = m[m._merge == "right_only"]["genome"].tolist()
if only_asm:
    print(f"\nWARNING: {len(only_asm)} assemblies have no fastbaps group, e.g. {only_asm[:5]}")
if only_fb:
    print(f"WARNING: {len(only_fb)} fastbaps entries have no assembly, e.g. {only_fb[:5]}")
if only_asm or only_fb:
    print("This is almost always an ID-formatting mismatch, not missing data. "
          "Fix it before continuing: unmatched genomes are silently dropped "
          "from every downstream model.\n")

m = m[m._merge == "both"].drop(columns="_merge")
m[["genome", "fastbaps_group", "n_contigs", "N50", "total_len"]].to_csv(
    out, sep="\t", index=False)

print(f"wrote {out}: {len(m)} genomes")
print(m.groupby("fastbaps_group").size().sort_values(ascending=False))
PY

rm -f .seqkit_stats.tsv
echo
echo "next: check the group labels above, then set --focal to the one you care about"
