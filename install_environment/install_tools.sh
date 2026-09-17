#!/usr/bin/env bash
# ==================================================================
# install_tools.sh
#
# Installs everything needed for the IS / transposition analyses.
# Idempotent: safe to re-run. Skips environments that already exist
# unless you pass --force, and skips databases already downloaded.
#
# Usage:
#   ./install_tools.sh                    # core env only
#   ./install_tools.sh --optional         # core + optional env
#   ./install_tools.sh --optional --dbs   # also download databases (~15 GB)
#   ./install_tools.sh --dbdir /scratch/username/dbs --optional --dbs
#   ./install_tools.sh --force            # rebuild environments
#   ./install_tools.sh --verify-only      # just check what works
#
#    chmod +x install_tools.sh
#   ./install_tools.sh --optional --dbs --dbdir /path/with/space
#
# On a cluster, run this on a compute node, not the login node. The
# solve is CPU and memory hungry and login-node limits will kill it:
#   srun -c 8 --mem 16G -t 4:00:00 --pty bash install_tools.sh --optional --dbs
# ==================================================================
set -uo pipefail

# ---------- configuration ----------------------------------------
CORE_ENV=ismobile
OPT_ENV=ismobile_extra
DBDIR="${PWD}/databases"
DO_OPTIONAL=0
DO_DBS=0
FORCE=0
VERIFY_ONLY=0
LOG="install_$(date +%Y%m%d_%H%M%S).log"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --optional)     DO_OPTIONAL=1; shift ;;
    --dbs)          DO_DBS=1; shift ;;
    --dbdir)        DBDIR="$2"; shift 2 ;;
    --force)        FORCE=1; shift ;;
    --verify-only)  VERIFY_ONLY=1; shift ;;
    -h|--help)      sed -n '2,25p' "$0"; exit 0 ;;
    *) echo "unknown option: $1"; exit 1 ;;
  esac
done

exec > >(tee -a "$LOG") 2>&1
say()  { printf '\n\033[1m==> %s\033[0m\n' "$*"; }
warn() { printf '\033[33m[warn] %s\033[0m\n' "$*"; }
fail() { printf '\033[31m[fail] %s\033[0m\n' "$*"; }
ok()   { printf '\033[32m[ ok ] %s\033[0m\n' "$*"; }

# ---------- locate conda -----------------------------------------
say "Locating conda"
if ! command -v conda >/dev/null 2>&1; then
  # common cluster pattern: conda lives behind a module
  if command -v module >/dev/null 2>&1; then
    module load miniconda3 2>/dev/null || module load anaconda3 2>/dev/null || true
  fi
fi
if ! command -v conda >/dev/null 2>&1; then
  fail "conda not found. Install miniforge first:"
  cat <<'EOF'
  curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
  bash Miniforge3-$(uname)-$(uname -m).sh -b -p "$HOME/miniforge3"
  "$HOME/miniforge3/bin/conda" init bash && exec bash
EOF
  exit 1
fi
CONDA_BASE=$(conda info --base)
source "${CONDA_BASE}/etc/profile.d/conda.sh"
ok "conda at ${CONDA_BASE}"

# mamba is not optional in practice. The optional env will not solve in
# reasonable time with the classic conda solver.
if command -v mamba >/dev/null 2>&1; then
  SOLVER=mamba
elif conda config --show solver 2>/dev/null | grep -q libmamba; then
  SOLVER=conda
  ok "using conda with the libmamba solver"
else
  say "Installing mamba into base (needed for a tractable solve)"
  conda install -y -n base -c conda-forge mamba || {
    warn "mamba install failed; falling back to conda + libmamba solver"
    conda install -y -n base -c conda-forge conda-libmamba-solver
    conda config --set solver libmamba
  }
  SOLVER=$(command -v mamba >/dev/null 2>&1 && echo mamba || echo conda)
fi
ok "solver: ${SOLVER}"

# ---------- environment creation ---------------------------------
env_exists() { conda env list | awk '{print $1}' | grep -qx "$1"; }

make_env () {
  local name="$1" yml="$2"
  if [[ ! -f "$yml" ]]; then
    fail "missing $yml (expected next to this script)"; return 1
  fi
  if env_exists "$name"; then
    if [[ $FORCE -eq 1 ]]; then
      say "Removing existing env ${name} (--force)"
      conda env remove -y -n "$name"
    else
      ok "env ${name} already exists, skipping (use --force to rebuild)"
      return 0
    fi
  fi
  say "Creating env ${name} from ${yml} (this takes a while)"
  $SOLVER env create -f "$yml" || { fail "could not create ${name}"; return 1; }
  ok "created ${name}"
}

if [[ $VERIFY_ONLY -eq 0 ]]; then
  make_env "$CORE_ENV" env_core.yml
  [[ $DO_OPTIONAL -eq 1 ]] && make_env "$OPT_ENV" env_optional.yml
fi

# ---------- databases --------------------------------------------
# Each of these is a separate multi-GB download that the conda package
# deliberately does not bundle. Set --dbdir to somewhere with space and
# no purge policy; scratch is usually wrong for these.
if [[ $DO_DBS -eq 1 && $VERIFY_ONLY -eq 0 ]]; then
  mkdir -p "$DBDIR"
  conda activate "$OPT_ENV" 2>/dev/null || { fail "optional env missing; rerun with --optional"; exit 1; }

  if [[ ! -d "${DBDIR}/genomad_db" ]]; then
    say "geNomad database (~1.5 GB)"
    genomad download-database "$DBDIR" && ok "genomad db" || warn "genomad db failed"
  else ok "genomad db present"; fi

  if [[ ! -d "${DBDIR}/checkv-db-v1.5" && -z "$(ls -d ${DBDIR}/checkv-db-* 2>/dev/null)" ]]; then
    say "CheckV database (~7 GB)"
    checkv download_database "$DBDIR" && ok "checkv db" || warn "checkv db failed"
  else ok "checkv db present"; fi

  say "DefenseFinder models"
  defense-finder update --models-dir "${DBDIR}/defense-finder" && ok "defense-finder models" \
    || warn "defense-finder update failed"

  say "PADLOC database"
  padloc --db-update && ok "padloc db" || warn "padloc db update failed"

  say "VirSorter2 database (~11 GB, slowest of the set)"
  if [[ ! -d "${DBDIR}/virsorter2_db" ]]; then
    virsorter setup -d "${DBDIR}/virsorter2_db" -j 4 && ok "virsorter2 db" \
      || warn "virsorter2 db failed (rerun the same command to resume)"
  else ok "virsorter2 db present"; fi

  cat > db_paths.env <<EOF
# source this before running the optional analyses
export GENOMAD_DB="${DBDIR}/genomad_db"
export CHECKVDB="\$(ls -d ${DBDIR}/checkv-db-* 2>/dev/null | head -1)"
export DEFENSE_FINDER_MODELS="${DBDIR}/defense-finder"
export VIRSORTER_DB="${DBDIR}/virsorter2_db"
EOF
  ok "wrote db_paths.env"
  conda deactivate
fi

# ---------- verification -----------------------------------------
PASS=0; MISS=0
check_cmd () {
  if command -v "$1" >/dev/null 2>&1; then ok "$1"; PASS=$((PASS+1))
  else fail "$1 not found"; MISS=$((MISS+1)); fi
}
check_py () {
  if python -c "import $1" 2>/dev/null; then ok "python: $1"; PASS=$((PASS+1))
  else fail "python module $1 missing"; MISS=$((MISS+1)); fi
}
check_r () {
  if Rscript -e "suppressMessages(library($1))" >/dev/null 2>&1; then
    ok "R: $1"; PASS=$((PASS+1))
  else fail "R package $1 missing"; MISS=$((MISS+1)); fi
}

say "Verifying core environment"
conda activate "$CORE_ENV" || { fail "cannot activate ${CORE_ENV}"; exit 1; }
for c in isescan.py blastn makeblastdb seqkit hmmsearch Rscript; do check_cmd "$c"; done
for m in pandas numpy scipy Bio ete3; do check_py "$m"; done
for p in data.table MASS ape phytools ggplot2; do check_r "$p"; done

# ISEScan needs FragGeneScan and its training data on PATH; this is the
# single most common install failure and it only shows up at runtime.
say "Smoke-testing ISEScan on a toy contig"
TMP=$(mktemp -d)
python - "$TMP" <<'PY'
import random, sys, os
random.seed(0)
seq = "".join(random.choice("ACGT") for _ in range(20000))
with open(os.path.join(sys.argv[1], "toy.fna"), "w") as f:
    f.write(">toy\n" + seq + "\n")
PY
if isescan.py --seqfile "$TMP/toy.fna" --output "$TMP/out" --nthread 1 >/dev/null 2>&1; then
  ok "ISEScan runs end to end"
else
  fail "ISEScan failed at runtime (usually FragGeneScan train files); see $LOG"
  warn "workaround: docker pull quay.io/biocontainers/isescan"
  MISS=$((MISS+1))
fi
rm -rf "$TMP"
conda deactivate

if [[ $DO_OPTIONAL -eq 1 ]] && env_exists "$OPT_ENV"; then
  say "Verifying optional environment"
  conda activate "$OPT_ENV"
  for c in digIS.py panISa.py ismap bwa-mem2 samtools minimap2 \
           genomad mob_recon virsorter checkv defense-finder padloc flye medaka; do
    check_cmd "$c"
  done
  conda deactivate
fi

# ---------- summary ----------------------------------------------
say "Summary"
echo "passed: ${PASS}   missing: ${MISS}"
echo "log:    ${LOG}"
cat <<EOF

Next steps:
  conda activate ${CORE_ENV}
  bash 01_annotate_IS.sh assemblies is_out 16
  Rscript 02_IS_load_and_recency.R
  python 03_insertion_site_polymorphism.py --asm assemblies --meta metadata.tsv
  python 04_IS_selection_and_recombination.py --focal G4
EOF
[[ $MISS -eq 0 ]] || warn "some tools are missing; the core scripts need the core env items only"
exit 0
