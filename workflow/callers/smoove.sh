#!/bin/bash
# smoove.sh: smoove 0.2.8 (LUMPY + SVTyper; brentp/smoove container), default settings, no exclude file.
#   smoove.sh call SPECIES SAMPLE                one sample: smoove call --genotype -> background "none"
#   smoove.sh merge SPECIES BACKGROUND SAMPLE|all  smoove's population workflow for one trio member against one
#       background (members as in delly2_merge.sh): smoove merge of the members' calls -> smoove genotype of
#       every member at the merged sites -> smoove paste -> the trio member's column
# -> $CALLS/SPECIES/smoove/BACKGROUND/SAMPLE/SAMPLE.vcf.gz; sample columns are named after the sample (the trio
# BAMs' SM tags are lane names). PRPOS/PREND are kept in the single-sample calls: smoove merge needs them.
# BREAK_REF (environment) replaces the species reference (depth series: chr7.fa).
# Env: delly2 (bcftools); smoove from $SMOOVE_SIF.
source "$(dirname "$0")/common.sh"
SIF=${SMOOVE_SIF:-$CNVB_DATA/tools/smoove_0.2.8.sif}
MODE=$1; species_settings "$2"; REF=${BREAK_REF:-$REF}
SM="apptainer exec -B $CNVB_DATA,$CNVB_RESULTS,$(dirname "$REF") $SIF smoove"
rename() { echo "$1" > "$3.name"; bcftools reheader -s "$3.name" -o "$3" "$2"; rm "$3.name"; }   # rename SAMPLE IN OUT
case $MODE in
call)
    S=$3; W=$WORK/smoove/$SPECIES/call/$S; O=$CALLS/$SPECIES/smoove/none/$S; mkdir -p "$W" "$O"
    $SM call --genotype -p "$THREADS" --name "$S" --fasta "$REF" --outdir "$W" "$CNVB_RESULTS/bam_bwa/$SPECIES/$S.bam"
    rename "$S" "$W/$S-smoove.genotyped.vcf.gz" "$O/$S.vcf.gz" && bcftools index -f "$O/$S.vcf.gz"
    rm -r "$W" ;;
merge)
    BG=$3; S=$4; W=$WORK/smoove/$SPECIES/$BG/$S; mkdir -p "$W" && cd "$W"
    if [ "$S" = all ]; then members "$SPECIES" "$BG" all > members.txt
    else { members "$SPECIES" "$BG" "$(query_group "$SPECIES" "$BG" "$S")"; echo "$S"; } > members.txt; fi
    $SM merge --name merged --fasta "$REF" --outdir . $(sed "s#.*#$CALLS/$SPECIES/smoove/none/&/&.vcf.gz#" members.txt)
    while read -r m; do
        $SM genotype --removepr -p 1 --name "$m" --outdir . --fasta "$REF" --vcf merged.sites.vcf.gz \
            "$CNVB_RESULTS/bam_bwa/$SPECIES/$m.bam"
        rename "$m" "$m-smoove.genotyped.vcf.gz" "rh_$m.vcf.gz" && bcftools index -t "rh_$m.vcf.gz"   # paste needs the index
    done < members.txt
    $SM paste --name cohort --outdir . $(sed 's#.*#rh_&.vcf.gz#' members.txt)
    for m in $([ "$S" = all ] && cat members.txt || echo "$S"); do
        O=$CALLS/$SPECIES/smoove/$BG/$m; mkdir -p "$O"
        bcftools view -s "$m" cohort.smoove.square.vcf.gz -Oz -o "$O/$m.vcf.gz"
    done
    cd .. && rm -r "$W" ;;
*) echo "usage: smoove.sh call|merge ..." >&2; exit 1 ;;
esac
