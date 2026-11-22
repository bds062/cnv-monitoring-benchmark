#!/bin/bash
# Human 1000 Genomes trios: CRAMs -> BAMs -> per-sample outputs of the five callers in $CALLS/human/.
#   run_all.sh [STAGE ...]     default: every stage, in order
# Each stage needs the previous one finished; with SUBMIT set, submit one stage at a time.
source "$(dirname "$0")/../callers/common.sh"
H=$(dirname "$0"); C=$H/../callers
species_settings human
ALL=$(samples human)

stage() {
    case $1 in
    reference) job preprocess bash "$H/prepare_reference.sh" ;;
    bam)       for s in $ALL; do job preprocess bash "$H/cram_to_bam.sh" "$s"; done ;;
    prepare)   for c in gatk cnvcaller cnvkit; do job $c bash "$C/$c.sh" prepare human; done ;;
    sample)    for s in $ALL; do
                   job delly2 bash "$C/delly2.sh" human "$s"
                   job delly2 bash "$C/smoove.sh" call human "$s"
                   job gatk bash "$C/gatk.sh" counts human "$s"
                   job cnvcaller bash "$C/cnvcaller.sh" rd human "$s"
                   job cnvkit bash "$C/cnvkit.sh" coverage human "$s"
               done ;;
    ploidy)    for g in $(groups human crosstrio9); do job gatk bash "$C/gatk.sh" ploidy human crosstrio9 $g; done ;;
    model)     for g in $(groups human crosstrio9); do for i in $(seq -f %02g 0 $((GATK_SHARDS - 1))); do
                   job gatk bash "$C/gatk.sh" model human crosstrio9 $g $i; done; done ;;
    call)      for t in $(groups human crosstrio9); do       # GATK-gCNV: each held-out trio in one CASE run
                   job gatk bash "$C/gatk.sh" call human crosstrio9 $(awk -F'\t' -v t="$t" '$1 == "human" && $2 == t {print $4}' "$SAMPLES/trios.tsv")
               done
               for s in $ALL; do
                   job cnvcaller bash "$C/cnvcaller.sh" call human crosstrio9 "$s"
                   job cnvkit bash "$C/cnvkit.sh" call human crosstrio3 "$s"
                   job delly2 bash "$C/delly2_merge.sh" human crosstrio9 "$s"
                   job delly2 bash "$C/smoove.sh" merge human crosstrio9 "$s"
               done ;;
    call_y200) for s in $ALL; do job delly2 env MERGE_Y=200 bash "$C/delly2_merge.sh" human crosstrio9 "$s"; done ;;   # -> crosstrio9_y200
    *) echo "unknown stage $1" >&2; exit 1 ;;
    esac
}
for st in ${@:-reference bam prepare sample ploidy model call call_y200}; do stage "$st"; done
