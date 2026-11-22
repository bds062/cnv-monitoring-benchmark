#!/bin/bash
# Time one benchmark step at the fixed configuration of the paper's run-time table and append one record to
# data/runtime/fixed_cores.tsv. Run it on a node with THREADS cores (8 in the paper), once per repeat:
#   bash workflow/runtime/timed.sh SPECIES STEP LABEL REP -- command [args...]
# The command is the step's line from the caller scripts in workflow/callers/ with -t/-p/--nproc set to $THREADS,
# calling the Trio1 child (zea) or NA19240 (human) against crosstrio-9 (CNVkit crosstrio-3), with earlier steps'
# outputs as inputs so that each step is timed on its own.
source "$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)/workflow/lib.sh"
SP=$1 STEP=$2 LABEL=$3 REP=$4; shift 4; [ "${1:-}" = "--" ] && shift
OUT=$CNVB_REPO/data/runtime/fixed_cores.tsv
U=$(mktemp); export OMP_NUM_THREADS=$THREADS
set +e; /usr/bin/time -v -o "$U" "$@"; RC=$?; set -e
CPU=$(lscpu | awk -F: '/Model name/ {gsub(/^ +/, "", $2); print $2; exit}')
awk -F': ' -v OFS='\t' -v sp="$SP" -v st="$STEP" -v lb="$LABEL" -v rp="$REP" -v cpu="$CPU" -v n="$THREADS" -v rc="$RC" \
    -v rs="${SLURM_RESTART_COUNT:-0}" '/Elapsed \(wall clock\)/ {w = $2} /User time/ {u = $2} /System time/ {s = $2}
    /Maximum resident set size/ {m = $2} END {print sp, st, lb, rp, cpu, n, rs, rc, w, u + s, m}' "$U" >> "$OUT"
rm -f "$U"; exit $RC
