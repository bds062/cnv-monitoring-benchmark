#!/bin/bash
# Chromosome 7 depth series (Additional file 1, Section S12), run stage by stage; each stage needs the earlier
# ones finished. Set SUBMIT to send every job to a scheduler (see workflow/lib.sh).
#   bash workflow/downsample/run_all.sh prepare   chromosome 7 bins and CNVcaller windows
#   bash workflow/downsample/run_all.sh bams      chromosome 7 BAMs of the 12 trio moths at full depth, 10x and 5x
#   bash workflow/downsample/run_all.sh moth      per moth and level: Delly2 discovery and delly cnv, smoove call,
#                                                 GATK-gCNV read counts, CNVcaller read depth, CNVkit coverage
#   bash workflow/downsample/run_all.sh models    GATK-gCNV crosstrio-9 models, per level and query trio
#   bash workflow/downsample/run_all.sh query     every moth against its background at its level: Delly2 and
#                                                 smoove merges, GATK-gCNV CASE, CNVcaller, CNVkit
#   bash workflow/downsample/run_all.sh pop       Delly2 copy number: delly cnv -l per moth, then the population
#   bash workflow/downsample/run_all.sh pop_call   read-depth workflow against crosstrio-9 (delly2_pop.sh)
# Needs workflow/zea through its 'sample' stage (Bowtie2 BAMs, GATK-gCNV bins, CNVcaller windows, CNVkit bins)
# and workflow/monitoring 'prepare' and 'chr7' (the trio moths' chromosome 7 BWA-MEM realignments).
# Delly2 and smoove run on the chromosome 7 BWA-MEM BAMs against the chromosome 7 reference (BREAK_REF=chr7.fa),
# and the Delly2 merge keeps every site (delly merge -y 0, as in every Delly2 run of the paper).
source "$(dirname "$0")/common.sh"
D=$(cd "$(dirname "$0")" && pwd)
TRIO=$(samples zea | tr '\n' ' ')
at() { echo "env CNVB_RESULTS=$DS/$1"; }       # at LEVEL: run a workflow/callers script on that level's tree
brk() { echo "$(at "$1") BREAK_REF=$REF7"; }    # brk LEVEL: the same, breakpoint callers on the chr7 reference

case ${1:-} in
prepare)
    job preprocess bash "$D/chr7.sh" prepare ;;
bams)
    for S in $TRIO; do job delly2 bash "$D/chr7.sh" bams "$S"; done ;;
moth)
    for L in "${LEVELS[@]}"; do for S in $TRIO; do
        job delly2 $(brk "$L") bash "$CALLERS/delly2.sh" zea "$S"
        job delly2 bash "$D/chr7.sh" delly_cnv "$L" "$S"
        job delly2 $(brk "$L") bash "$CALLERS/smoove.sh" call zea "$S"
        job gatk bash "$D/chr7.sh" gatk_counts "$L" "$S"
        job cnvcaller $(at "$L") bash "$CALLERS/cnvcaller.sh" rd zea "$S"
        job cnvkit bash "$D/chr7.sh" cnvkit_coverage "$L" "$S"
    done; done ;;
models)
    for L in "${LEVELS[@]}"; do for T in Trio1 Trio2 Trio3 Trio4; do job gatk bash "$D/chr7.sh" gatk_model "$L" "$T"; done; done ;;
query)
    for L in "${LEVELS[@]}"; do for S in $TRIO; do
        job delly2 $(brk "$L") MERGE_Y=0 bash "$CALLERS/delly2_merge.sh" zea crosstrio9 "$S"
        job delly2 $(brk "$L") bash "$CALLERS/smoove.sh" merge zea crosstrio9 "$S"
        job gatk bash "$D/chr7.sh" gatk_case "$L" "$S"
        job cnvcaller $(at "$L") bash "$CALLERS/cnvcaller.sh" call zea crosstrio9 "$S"
        job cnvkit bash "$D/chr7.sh" cnvkit_call "$L" "$S"
    done; done ;;
pop)
    for L in "${LEVELS[@]}"; do for S in $TRIO; do
        job delly2 $(brk "$L") MAP=$MAP7 bash "$CALLERS/delly2_pop.sh" cnvl zea "$S"
    done; done ;;
pop_call)
    for L in "${LEVELS[@]}"; do for S in $TRIO; do
        job delly2 $(brk "$L") MAP=$MAP7 bash "$CALLERS/delly2_pop.sh" panel zea crosstrio9 "$S"
    done; done ;;
*)
    sed -n '2,17p' "$0"; exit 1 ;;
esac
