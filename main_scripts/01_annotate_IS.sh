#!/usr/bin/env bash
# ------------------------------------------------------------------
# 01_annotate_IS.sh
# Proper IS-element annotation (not just "putative transposase" from
# Prokka/Bakta), extraction of IS copies, and all-vs-all identity so
# you can date transposition bursts.
#
# Inputs:
#   ASM_DIR   directory of assemblies, one per genome, <genome>.fna
#   META      TSV: genome, fastbaps_group, n_contigs, N50, total_len
# Outputs:
#   isescan_out/          raw ISEScan output
#   IS_calls.tsv          tidy per-copy table for all genomes
#   IS_seqs.fna           nucleotide sequence of every IS copy
#   IS_identity.tsv       all-vs-all % identity within IS families
#
# conda create -n ismobile -c bioconda isescan blast seqkit diamond
# ------------------------------------------------------------------
set -euo pipefail

ASM_DIR=${1:-assemblies}
OUT=${2:-is_out}
THREADS=${3:-16}

mkdir -p "$OUT"/isescan_out

# ---- 1. ISEScan on every assembly -------------------------------
# ISEScan finds full IS elements: transposase ORF + terminal inverted
# repeats (TIR) + target site duplication (TSD). The TIR/TSD columns
# are what let you separate intact, recently active copies from
# degenerate relics, which a Prokka "transposase" label cannot do.
ls "$ASM_DIR"/*.fna | xargs -P "$THREADS" -I{} bash -c '
  g=$(basename {} .fna)
  isescan.py --seqfile {} --output '"$OUT"'/isescan_out/$g --nthread 1
'

# ---- 2. Collate into one tidy table ------------------------------
python3 - "$OUT" <<'PY'
import sys, glob, os, pandas as pd
out = sys.argv[1]
rows = []
for f in glob.glob(os.path.join(out, "isescan_out", "*", "*.raw")) + \
         glob.glob(os.path.join(out, "isescan_out", "*", "*.tsv")):
    if f.endswith(".tsv") is False:
        continue
    g = os.path.basename(os.path.dirname(f))
    try:
        d = pd.read_csv(f, sep="\t")
    except Exception:
        continue
    if d.empty:
        continue
    d["genome"] = g
    rows.append(d)
tab = pd.concat(rows, ignore_index=True)

# ISEScan column names vary slightly by version; normalise the ones we need
ren = {"seqID": "contig", "isBegin": "start", "isEnd": "end",
       "isLen": "is_len", "type": "is_type", "irId": "ir_id",
       "irLen": "ir_len", "orfLen": "orf_len"}
tab = tab.rename(columns={k: v for k, v in ren.items() if k in tab.columns})

# complete = ISEScan type "c": transposase ORF plus both TIRs present
tab["complete"] = tab.get("is_type", "").astype(str).str.startswith("c")
# a plausible-length ORF is the second filter for "still functional"
tab["intact_orf"] = tab["orf_len"] > 0.8 * tab.groupby("family")["orf_len"].transform("max")

tab.to_csv(os.path.join(out, "IS_calls.tsv"), sep="\t", index=False)
print(tab.groupby(["family", "complete"]).size())
PY

# ---- 3. Pull the nucleotide sequence of every copy ---------------
python3 - "$OUT" "$ASM_DIR" <<'PY'
import sys, os, pandas as pd
from Bio import SeqIO
out, asm = sys.argv[1], sys.argv[2]
tab = pd.read_csv(os.path.join(out, "IS_calls.tsv"), sep="\t")
seqs = {}
with open(os.path.join(out, "IS_seqs.fna"), "w") as fh:
    for g, sub in tab.groupby("genome"):
        idx = SeqIO.index(os.path.join(asm, f"{g}.fna"), "fasta")
        for i, r in sub.iterrows():
            s = idx[r["contig"]].seq[int(r["start"]) - 1:int(r["end"])]
            name = f"{g}|{r['contig']}|{int(r['start'])}|{int(r['end'])}|{r['family']}"
            fh.write(f">{name}\n{s}\n")
PY

# ---- 4. All-vs-all identity within families ----------------------
# Recent transposition bursts show up as a pile-up of near-identical
# copies. Old, vertically inherited copies diverge with the core
# genome, so their identity distribution tracks core divergence.
# Compare that distribution in the focal group vs the other groups.
makeblastdb -in "$OUT"/IS_seqs.fna -dbtype nucl -out "$OUT"/IS_db >/dev/null
blastn -query "$OUT"/IS_seqs.fna -db "$OUT"/IS_db \
       -outfmt '6 qseqid sseqid pident length qlen slen' \
       -perc_identity 70 -num_threads "$THREADS" \
       -max_target_seqs 100000 > "$OUT"/IS_identity_raw.tsv

awk -F'\t' 'BEGIN{OFS="\t"; print "q","s","pident","alen","qlen","slen"}
     $1!=$2 && $4 > 0.8*$5 {print}' "$OUT"/IS_identity_raw.tsv > "$OUT"/IS_identity.tsv

echo "done: $OUT/IS_calls.tsv, $OUT/IS_seqs.fna, $OUT/IS_identity.tsv"
