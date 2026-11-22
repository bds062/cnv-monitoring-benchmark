#!/bin/bash
# Make every table and figure of the paper from the caller outputs in $CNVB_RESULTS, in dependency order.
#   bash analysis/make_all.sh            (from the repository root, in the envs/analysis.yml environment)
# Tables go to $CNVB_RESULTS/tables/, figures to $CNVB_RESULTS/figures/. README.md maps each paper figure and
# table to its script and output file. Stops at the first script that fails.
set -euo pipefail
cd "$(dirname "$0")"
export MPLBACKEND=Agg
for s in trio_benchmark unrelated_parents blacklist_overlap ddpcr human_truthset data_quality coverage_plots \
         tryp77_region delly2_extras cnvkit_depth reference_stats runtime datasets chromosome_map \
         monitoring_panels monitoring_query monitoring_scan monitoring_sweep cnvcaller_f0 downsample_chr7 \
         overview_figures; do
    echo "== $s"
    python "$s.py"
done
