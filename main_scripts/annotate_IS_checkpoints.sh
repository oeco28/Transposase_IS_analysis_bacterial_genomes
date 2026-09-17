#!/usr/bin/env bash
# ==================================================================
# 01_annotate_IS.sh  (v2 - resumable)
#
# v1 restarted everything on every invocation. At ~1900 genomes a
# wall-clock timeout meant losing the whole run. This version
# checkpoints at four levels:
#
#   1. per genome   - a genome with a valid .done marker is skipped
#   2. per stage    - annotate / collate / extract / blast run only if
#                     their output is missing or older than their input
#   3. per BLAST chunk - the all-vs-all is split into chunks, each with
#                     its own marker, so a killed BLAST resumes
#   4. atomically   - markers are written only after the output is
#                     validated, so a job killed mid-genome is redone
#                     rather than silently left truncated
#
# Usage:
#   bash 01_annotate_IS.sh assemblies is_out 16
#   bash 01_annotate_IS.sh assemblies is_out 16 --force        # redo all
#   bash 01_annotate_IS.sh assemblies is_out 16 --stage blast  # one stage
#   bash 01_annotate_IS.sh assemblies is_out 16 --status       # what is left
#
# SLURM array (the right way to do 1900 genomes):
#   bash 01_annotate_IS.sh assemblies is_out 1 --status   # writes todo.txt
#   sbatch --array=1-$(wc -l < is_out/todo.txt)%50 array_job.sh
#   # then rerun this script normally to do the collate/extract/blast
#   # stages once the array has finished.
# ==================================================================
set -uo pipefail

ASM_DIR=${1:-assemblies}
OUT=${2:-is_out}
THREADS=${3:-16}
shift 3 2>/dev/null || true

FORCE=0; ONLY_STAGE=""; STATUS_ONLY=0
while [[ $# -gt 0 ]]; do
  case "$1" in
    --force)  FORCE=1; shift ;;
    --stage)  ONLY_STAGE="$2"; shift 2 ;;
    --status) STATUS_ONLY=1; shift ;;
    *) echo "unknown option: $1"; exit 1 ;;
  esac
done

RAW="$OUT/isescan_out"
MARK="$OUT/.markers"
mkdir -p "$RAW" "$MARK"

run_stage () { [[ -z "$ONLY_STAGE" || "$ONLY_STAGE" == "$1" ]]; }
newer_than () { [[ "$1" -nt "$2" ]]; }

# ---- inventory ---------------------------------------------------
mapfile -t ALL < <(ls "$ASM_DIR"/*.fna 2>/dev/null | sort)
[[ ${#ALL[@]} -gt 0 ]] || { echo "no *.fna in $ASM_DIR"; exit 1; }

TODO=(); DONE=0
for f in "${ALL[@]}"; do
  g=$(basename "$f" .fna)
  if [[ $FORCE -eq 0 && -f "$MARK/$g.done" ]]; then
    DONE=$((DONE+1))
  else
    TODO+=("$f")
  fi
done

echo "genomes total: ${#ALL[@]}   already done: ${DONE}   remaining: ${#TODO[@]}"

if [[ $STATUS_ONLY -eq 1 ]]; then
  printf '%s\n' "${TODO[@]}" > "$OUT/todo.txt"
  echo "wrote $OUT/todo.txt (${#TODO[@]} genomes) for use as a SLURM array input"
  # emit a ready-to-submit array script alongside it
  cat > array_job.sh <<EOF
#!/usr/bin/env bash
#SBATCH --job-name=isescan
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=2:00:00
#SBATCH --output=logs/isescan_%A_%a.out
set -euo pipefail
mkdir -p logs
source "\$(conda info --base)/etc/profile.d/conda.sh"
conda activate ismobile
f=\$(sed -n "\${SLURM_ARRAY_TASK_ID}p" $OUT/todo.txt)
g=\$(basename "\$f" .fna)
tmp=\$(mktemp -d)
if isescan.py --seqfile "\$f" --output "\$tmp/\$g" --nthread 2; then
  rm -rf "$RAW/\$g"
  mv "\$tmp/\$g" "$RAW/\$g"
  touch "$MARK/\$g.done"
else
  echo "FAILED: \$g" >&2; exit 1
fi
rm -rf "\$tmp"
EOF
  chmod +x array_job.sh
  echo "wrote array_job.sh"
  exit 0
fi

# ---- stage 1: ISEScan, per-genome checkpointed -------------------
# Each genome runs into a temp dir and is only moved into place and
# marked done after ISEScan exits cleanly. A job killed mid-genome
# leaves no marker, so the next run redoes that genome and only that
# genome. Without this you get a truncated .tsv that looks valid.
if run_stage annotate && [[ ${#TODO[@]} -gt 0 ]]; then
  echo "=== stage 1: ISEScan on ${#TODO[@]} genomes ==="
  printf '%s\n' "${TODO[@]}" | xargs -P "$THREADS" -I{} bash -c '
    f="{}"; g=$(basename "$f" .fna)
    tmp=$(mktemp -d)
    if isescan.py --seqfile "$f" --output "$tmp/$g" --nthread 1 >/dev/null 2>&1; then
      rm -rf "'"$RAW"'/$g"
      mv "$tmp/$g" "'"$RAW"'/$g"
      touch "'"$MARK"'/$g.done"
    else
      echo "FAILED: $g" >> "'"$OUT"'/failed_genomes.txt"
    fi
    rm -rf "$tmp"
  '
  if [[ -s "$OUT/failed_genomes.txt" ]]; then
    n=$(sort -u "$OUT/failed_genomes.txt" | wc -l)
    echo "WARNING: $n genomes failed; see $OUT/failed_genomes.txt"
    echo "  rerun this script to retry them (successful ones are skipped)"
  fi
elif run_stage annotate; then
  echo "=== stage 1: all genomes already annotated, skipping ==="
fi

# ---- stage 2: collate -------------------------------------------
# Cheap, but redone whenever a new genome finished, so gate it on the
# newest marker rather than on existence alone.
CALLS="$OUT/IS_calls.tsv"
NEWEST_MARK=$(ls -t "$MARK"/*.done 2>/dev/null | head -1)
if run_stage collate && { [[ $FORCE -eq 1 ]] || [[ ! -f "$CALLS" ]] || \
   { [[ -n "$NEWEST_MARK" ]] && newer_than "$NEWEST_MARK" "$CALLS"; }; }; then
  echo "=== stage 2: collating IS calls ==="
  python3 - "$OUT" <<'PY'
import sys, glob, os
import pandas as pd

out = sys.argv[1]
rows = []
for f in glob.glob(os.path.join(out, "isescan_out", "*", "**", "*.tsv"),
                   recursive=True):
    g = os.path.basename(os.path.dirname(f))
    # ISEScan nests output under the input path; walk up to the genome dir
    parts = os.path.relpath(f, os.path.join(out, "isescan_out")).split(os.sep)
    g = parts[0]
    try:
        d = pd.read_csv(f, sep="\t")
    except Exception:
        continue
    if d.empty or "seqID" not in d.columns and "family" not in d.columns:
        continue
    d["genome"] = g
    rows.append(d)

if not rows:
    sys.exit("no ISEScan tables found; check is_out/isescan_out/")

tab = pd.concat(rows, ignore_index=True)
ren = {"seqID": "contig", "isBegin": "start", "isEnd": "end", "isLen": "is_len",
       "type": "is_type", "irId": "ir_id", "irLen": "ir_len", "orfLen": "orf_len"}
tab = tab.rename(columns={k: v for k, v in ren.items() if k in tab.columns})
tab["complete"] = tab.get("is_type", "").astype(str).str.startswith("c")
if "orf_len" in tab.columns and "family" in tab.columns:
    tab["intact_orf"] = tab["orf_len"] > 0.8 * tab.groupby("family")["orf_len"].transform("max")

tmp = os.path.join(out, ".IS_calls.tmp")
tab.to_csv(tmp, sep="\t", index=False)
os.replace(tmp, os.path.join(out, "IS_calls.tsv"))   # atomic
print(f"{len(tab)} IS copies across {tab.genome.nunique()} genomes")
print(tab.groupby(["family", "complete"]).size().head(20))
PY
elif run_stage collate; then
  echo "=== stage 2: IS_calls.tsv up to date, skipping ==="
fi

# ---- stage 3: extract sequences ---------------------------------
SEQS="$OUT/IS_seqs.fna"
if run_stage extract && { [[ $FORCE -eq 1 ]] || [[ ! -f "$SEQS" ]] || \
   newer_than "$CALLS" "$SEQS"; }; then
  echo "=== stage 3: extracting IS sequences ==="
  python3 - "$OUT" "$ASM_DIR" <<'PY'
import sys, os
import pandas as pd
from Bio import SeqIO

out, asm = sys.argv[1], sys.argv[2]
tab = pd.read_csv(os.path.join(out, "IS_calls.tsv"), sep="\t")
tmp = os.path.join(out, ".IS_seqs.tmp")
n = 0
with open(tmp, "w") as fh:
    for g, sub in tab.groupby("genome"):
        p = os.path.join(asm, f"{g}.fna")
        if not os.path.exists(p):
            continue
        idx = SeqIO.index(p, "fasta")
        for _, r in sub.iterrows():
            if r["contig"] not in idx:
                continue
            s = idx[r["contig"]].seq[int(r["start"]) - 1:int(r["end"])]
            fh.write(f">{g}|{r['contig']}|{int(r['start'])}|{int(r['end'])}|"
                     f"{r['family']}\n{s}\n")
            n += 1
os.replace(tmp, os.path.join(out, "IS_seqs.fna"))
print(f"wrote {n} IS sequences")
PY
elif run_stage extract; then
  echo "=== stage 3: IS_seqs.fna up to date, skipping ==="
fi

# ---- stage 4: chunked all-vs-all BLAST --------------------------
# At ~1900 genomes this is tens of thousands of sequences and is by far
# the longest step. Split the query into chunks so a killed job resumes
# at the chunk boundary instead of restarting the whole search.
IDENT="$OUT/IS_identity.tsv"
if run_stage blast && { [[ $FORCE -eq 1 ]] || [[ ! -f "$IDENT" ]] || \
   newer_than "$SEQS" "$IDENT"; }; then
  echo "=== stage 4: all-vs-all identity ==="
  CHUNKDIR="$OUT/.blast_chunks"
  mkdir -p "$CHUNKDIR"

  if [[ ! -f "$OUT/IS_db.nsq" ]] || newer_than "$SEQS" "$OUT/IS_db.nsq"; then
    makeblastdb -in "$SEQS" -dbtype nucl -out "$OUT/IS_db" >/dev/null
  fi

  if [[ $FORCE -eq 1 ]] || [[ ! -f "$CHUNKDIR/.split.done" ]] || \
     newer_than "$SEQS" "$CHUNKDIR/.split.done"; then
    rm -f "$CHUNKDIR"/chunk_* 2>/dev/null
    seqkit split2 -p 40 -O "$CHUNKDIR" -f "$SEQS" >/dev/null 2>&1 || \
      awk -v d="$CHUNKDIR" 'BEGIN{n=0;c=0}
           /^>/{if(n%2000==0)c++; n++} {print > (d"/chunk_"c".fna")}' "$SEQS"
    touch "$CHUNKDIR/.split.done"
  fi

  mapfile -t CHUNKS < <(ls "$CHUNKDIR"/*.fna 2>/dev/null | sort)
  echo "${#CHUNKS[@]} chunks"
  for c in "${CHUNKS[@]}"; do
    base=$(basename "$c" .fna)
    if [[ $FORCE -eq 0 && -f "$CHUNKDIR/$base.done" ]]; then continue; fi
    echo "  blast $base"
    if blastn -query "$c" -db "$OUT/IS_db" \
              -outfmt '6 qseqid sseqid pident length qlen slen' \
              -perc_identity 70 -num_threads "$THREADS" \
              -max_target_seqs 100000 > "$CHUNKDIR/$base.raw.tmp"; then
      mv "$CHUNKDIR/$base.raw.tmp" "$CHUNKDIR/$base.raw"
      touch "$CHUNKDIR/$base.done"
    else
      echo "  chunk $base failed; rerun to resume here"; exit 1
    fi
  done

  { echo -e "q\ts\tpident\talen\tqlen\tslen"
    cat "$CHUNKDIR"/*.raw | awk -F'\t' 'BEGIN{OFS="\t"} $1!=$2 && $4 > 0.8*$5'
  } > "$OUT/.IS_identity.tmp"
  mv "$OUT/.IS_identity.tmp" "$IDENT"
  echo "wrote $IDENT ($(($(wc -l < "$IDENT") - 1)) pairs)"
elif run_stage blast; then
  echo "=== stage 4: IS_identity.tsv up to date, skipping ==="
fi

echo
echo "outputs: $CALLS, $SEQS, $IDENT"
echo "rerun this script any time; completed work is skipped"
