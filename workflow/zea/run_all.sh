#!/bin/bash
# H. zea: public reads -> BAMs -> per-sample outputs of the five callers in $CALLS/zea/.
#   run_all.sh [STAGE ...]     default: every stage, in order
# Each stage needs the previous one finished; with SUBMIT set, submit one stage at a time.
source "$(dirname "$0")/../callers/common.sh"
Z=$(dirname "$0"); C=$Z/../callers
species_settings zea
TRIO=$(samples zea)
PUB=$(awk -F'\t' 'NR > 1 {print $1}' "$SAMPLES/published_moths.tsv")                  # 41 published moths
PUB_BG=$(awk -F'\t' 'NR == FNR {if (FNR > 1) p[$1]; next} $1 == "zea" && ($4 in p) {print $4}' \
    "$SAMPLES/published_moths.tsv" "$SAMPLES/backgrounds.tsv" | sort -u)                 # the 25 in backgrounds
TWICE=$(awk -F'\t' 'NR > 1 && ($4 == "LA2002" || $4 == "Benzon") {print $1}' "$SAMPLES/published_moths.tsv")
ALL="$TRIO $PUB"
BREAK="$TRIO $PUB_BG"           # moths of the breakpoint callers (BWA-MEM), as in the paper's runs
GATK_BG="crosstrio9 pooled12 susceptible13 neutral18 neutral9"
CNVCALLER_BG="crosstrio9 pooled12 susceptible13 neutral18 neutral9"
CNVKIT_BG="crosstrio3 crosstrio3_auto pooled12 susceptible13 neutral18"
BREAK_BG="crosstrio9 pooled12 susceptible13 neutral18 neutral9"       # Delly2 merge and smoove

stage() {
    case $1 in
    reference) job preprocess bash "$Z/prepare_reference.sh" ;;
    align)     for s in $ALL; do job preprocess bash "$Z/align.sh" "$s"; done
               for s in $TWICE; do job preprocess bash "$Z/align.sh" "$s" taylor2021; done ;;   # second BAM, lib.sh
    speedseq)  for s in $BREAK; do job speedseq bash "$Z/speedseq.sh" "$s"; done ;;
    prepare)   job delly2 bash "$C/delly2_cnv.sh" prepare zea
               for c in gatk cnvcaller cnvkit; do job $c bash "$C/$c.sh" prepare zea; done ;;
    sample)    for s in $TRIO; do job delly2 bash "$C/delly2_cnv.sh" call zea "$s"; done
               for s in $BREAK; do
                   job delly2 bash "$C/delly2.sh" zea "$s"
                   job delly2 bash "$C/smoove.sh" call zea "$s"
               done
               for s in $ALL; do        # all 47: the monitoring analyses reuse the counts and read depth
                   job gatk bash "$C/gatk.sh" counts zea "$s"
                   job cnvcaller bash "$C/cnvcaller.sh" rd zea "$s"
                   job cnvkit bash "$C/cnvkit.sh" coverage zea "$s"
               done ;;
    ploidy)    for b in $GATK_BG; do for g in $(groups zea $b); do job gatk bash "$C/gatk.sh" ploidy zea $b $g; done; done ;;
    model)     for b in $GATK_BG; do for g in $(groups zea $b); do for i in $(seq -f %02g 0 $((GATK_SHARDS - 1))); do
                   job gatk bash "$C/gatk.sh" model zea $b $g $i; done; done; done ;;
    call)      for s in $TRIO; do
                   for b in $GATK_BG; do job gatk bash "$C/gatk.sh" call zea $b "$s"; done
                   for b in $CNVCALLER_BG; do job cnvcaller bash "$C/cnvcaller.sh" call zea $b "$s"; done
                   for b in $CNVKIT_BG; do job cnvkit bash "$C/cnvkit.sh" call zea $b "$s"; done
               done
               for b in $BREAK_BG; do                    # pooled12 is one run over the 12 trio moths
                   for s in $([ $b = pooled12 ] && echo all || echo $TRIO); do
                       job delly2 bash "$C/delly2_merge.sh" zea $b "$s"
                       job delly2 bash "$C/smoove.sh" merge zea $b "$s"
                   done
               done ;;
    call_y200) for b in $BREAK_BG; do                    # Delly2 merge with its default minimum site quality -> BG_y200
                   for s in $([ $b = pooled12 ] && echo all || echo $TRIO); do
                       job delly2 env MERGE_Y=200 bash "$C/delly2_merge.sh" zea $b "$s"
                   done
               done ;;
    pop)       for s in $TRIO $PUB_BG; do job delly2 bash "$C/delly2_pop.sh" cnvl zea "$s"; done ;;
    pop_call)  for s in $TRIO; do job delly2 bash "$C/delly2_pop.sh" panel zea crosstrio9 "$s"; done   # genome-wide
               for b in pooled12 susceptible13 neutral18 neutral9; do                                # Tryp77 only
                   for s in $([ $b = pooled12 ] && echo all || echo $TRIO); do
                       job delly2 env POP_CHR7=1 bash "$C/delly2_pop.sh" panel zea $b "$s"
                   done
               done ;;
    *) echo "unknown stage $1" >&2; exit 1 ;;
    esac
}
for st in ${@:-reference align speedseq prepare sample ploidy model call call_y200 pop pop_call}; do stage "$st"; done
