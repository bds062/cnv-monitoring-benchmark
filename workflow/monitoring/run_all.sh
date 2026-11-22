#!/bin/bash
# Monitoring analyses on the published Louisiana time series, run stage by stage; each stage needs
# the earlier ones finished. Set SUBMIT to send every job to a scheduler (see workflow/lib.sh).
#   bash workflow/monitoring/run_all.sh prepare   read-depth truth of the 35 published moths, chromosome 7 reference,
#                                                 GATK-gCNV read counts and CNVcaller read depth of the 13 LA2002
#                                                 moths from their Taylor 2021-trimmed BAMs (in $T21, common.sh)
#   bash workflow/monitoring/run_all.sh chr7      chromosome 7 realignment and per-moth Delly2, 47 moths
#   bash workflow/monitoring/run_all.sh panels    time-series panels and the nested constant-size panels
#                                                 (ten series, 0-8 carriers) through the three callers;
#                                                 Delly2 also genotypes the 12 trio moths at the sites of the
#                                                 growing panels and of crosstrio-9
#   bash workflow/monitoring/run_all.sh query     the 12 trio moths against the growing panels (fixed query)
#                                                 and the nested constant-size panels (carrier titration)
#   bash workflow/monitoring/run_all.sh sweep     CNVcaller frequency sweep
#   bash workflow/monitoring/run_all.sh f0        every CNVcaller panel and query of the stages above again with
#                                                 -f 0 instead of 0.1, into $CNVB_RESULTS_F0/monitoring (default:
#                                                 results_f0 next to $CNVB_RESULTS; analysis/cnvcaller_f0.py)
# Needs workflow/zea run through its 'sample' stage (BAMs, GATK-gCNV counts and CNVcaller read depth of
# all 47 moths, the Delly2 mappability map) and, for crosstrio-9, its 'call' stage.
source "$(dirname "$0")/common.sh"
D=$(cd "$(dirname "$0")" && pwd)
TRIO=$(samples zea | tr '\n' ' ')
trio() { awk -F'\t' -v t="$1" '$1 == "zea" && $2 == t {printf "%s ", $4}' "$SAMPLES/trios.tsv"; }

case ${1:-} in
prepare)
    job delly2 bash "$D/read_depth.sh"
    job delly2 bash "$D/delly_chr7.sh" ref
    # $T21 is a results tree of its own: the main tree's bins and CNVcaller windows, its own counts and read depth
    mkdir -p "$T21/work/gatk/zea" "$T21/work/cnvcaller/zea"
    ln -sf "$GATK_ZEA/bins.interval_list" "$T21/work/gatk/zea/"
    ln -sf "$CNVC_ZEA/referenceDB.800" "$CNVC_ZEA/all_windows.link" "$T21/work/cnvcaller/zea/"
    for S in $(awk -F'\t' 'NR > 1 && $4 == "LA2002" {print $1}' "$SAMPLES/published_moths.tsv"); do
        job gatk env CNVB_RESULTS="$T21" bash "$D/../callers/gatk.sh" counts zea "$S"
        job cnvcaller env CNVB_RESULTS="$T21" bash "$D/../callers/cnvcaller.sh" rd zea "$S"
    done ;;
chr7)
    for S in $TRIO $(panel_members '2002+2012+2017'); do job delly2 bash "$D/delly_chr7.sh" "$S"; done ;;
panels)
    for P in $(panels "$SERIES|nest_.*"); do
        job gatk bash "$D/gatk_panel.sh" "$P"
        job cnvcaller bash "$D/cnvcaller_panel.sh" "$P"
        job delly2 bash "$D/delly_panel.sh" "$P" 200
        if [[ $P =~ ^($GROWING)$ ]]; then job delly2 bash "$D/delly_panel.sh" "$P" 0 $TRIO
        else job delly2 bash "$D/delly_panel.sh" "$P" 0; fi
    done
    for T in Trio1 Trio2 Trio3 Trio4; do job delly2 bash "$D/delly_panel.sh" "crosstrio9_$T" 0 $(trio $T); done ;;
query)
    for S in $TRIO; do
        for P in $(panels "$GROWING"); do
            job gatk bash "$D/gatk_panel.sh" "$P" "$S"
            job cnvcaller bash "$D/cnvcaller_panel.sh" "$P" "$S"
        done
        for P in $(panels 'nest_.*'); do
            job gatk bash "$D/gatk_panel.sh" "$P" "$S"         # one CASE run per job, see gatk_panel.sh
            job cnvcaller bash "$D/cnvcaller_panel.sh" "$P" "$S" chr7
        done
    done ;;
sweep)
    for P in $(panels 'sweep_.*'); do job cnvcaller bash "$D/cnvcaller_panel.sh" "$P"; done ;;
f0)
    F0="${CNVB_RESULTS_F0:-$(dirname "$CNVB_RESULTS")/results_f0}/monitoring"
    for P in $(panels "$SERIES|nest_.*|sweep_.*"); do job cnvcaller env CNVC_F=0 CNVC_OUT="$F0" bash "$D/cnvcaller_panel.sh" "$P"; done
    for S in $TRIO; do
        for P in $(panels "$GROWING"); do job cnvcaller env CNVC_F=0 CNVC_OUT="$F0" bash "$D/cnvcaller_panel.sh" "$P" "$S"; done
        for P in $(panels 'nest_.*'); do job cnvcaller env CNVC_F=0 CNVC_OUT="$F0" bash "$D/cnvcaller_panel.sh" "$P" "$S" chr7; done
    done ;;
*)
    sed -n '2,19p' "$0"; exit 1 ;;
esac
