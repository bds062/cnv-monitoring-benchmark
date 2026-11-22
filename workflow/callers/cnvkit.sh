#!/bin/bash
# CNVkit 0.9.13, WGS mode, default segmentation and call thresholds. Env: cnvkit.
#   cnvkit.sh prepare SPECIES                 bins: 5 kb over the accessible genome (flat batch, no samples);
#                                             H. zea also the 4,589 bp bins of the earlier protocol
#   cnvkit.sh coverage SPECIES SAMPLE         coverage of one sample on every bin set
#   cnvkit.sh call SPECIES BACKGROUND SAMPLE  reference from the background, then fix, segment, call, export
# crosstrio3 uses the protocol shared by both species (5 kb bins; human --diploid-parx-genome grch38 at
# every step). The other H. zea backgrounds (crosstrio3_auto = crosstrio3 members, pooled12,
# susceptible13, neutral18) used the earlier H. zea protocol, whose automatic bin size was 4,589 bp.
# In H. zea, chromosome 31 (Z) is named X.
source "$(dirname "$0")/common.sh"
species_settings "$2"
K=$WORK/cnvkit/$SPECIES
[ "$SPECIES" = zea ] && BINSETS="bins5k bins4589" || BINSETS=bins5k

case $1 in
prepare)
    for b in $BINSETS; do
        [ "$b" = bins4589 ] && SIZE=(--target-avg-size 4589) || SIZE=()
        mkdir -p "$K/$b" && cd "$K/$b"
        cnvkit.py batch --normal --fasta "$REF_DEPTH" --method wgs "${SIZE[@]}" "${PARX[@]}" \
            --output-reference flat_reference.cnn --output-dir . -p "$THREADS"
    done ;;
coverage)
    S=$3
    if [ "$SPECIES" = zea ]; then BAM=$(bam_x "$S" "$K"); else BAM=$CNVB_RESULTS/bam/human/$S.bam; fi
    for b in $BINSETS; do
        cnvkit.py coverage "$BAM" "$K/$b/"*.target.bed -p "$THREADS" -o "$K/$b/$S.targetcoverage.cnn"
        cnvkit.py coverage "$BAM" "$K/$b/"*.antitarget.bed -p 1 -o "$K/$b/$S.antitargetcoverage.cnn"
    done
    [ "$SPECIES" = human ] || rm "$BAM" "$BAM.bai" ;;
call)
    BG=$3; S=$4; O=$CALLS/$SPECIES/cnvkit/$BG/$S
    [ "$BG" = crosstrio3 ] && B=$K/bins5k || B=$K/bins4589
    mkdir -p "$O" && cd "$O"
    REFS=(); for m in $(members "$SPECIES" "$BG" "$(query_group "$SPECIES" "$BG" "$S")"); do
        REFS+=("$B/$m.targetcoverage.cnn" "$B/$m.antitargetcoverage.cnn"); done
    cnvkit.py reference "${REFS[@]}" -f "$REF_DEPTH" "${PARX[@]}" -o reference.cnn
    cnvkit.py fix "$B/$S.targetcoverage.cnn" "$B/$S.antitargetcoverage.cnn" reference.cnn "${PARX[@]}" -o "$S.cnr"
    cnvkit.py segment "$S.cnr" "${PARX[@]}" -p 4 -o "$S.cns"
    cnvkit.py call "$S.cns" "${PARX[@]}" -o "$S.call.cns"
    cnvkit.py export vcf "$S.call.cns" "${PARX[@]}" -i "$S" -o "$S.call.vcf"
    bgzip -f "$S.call.vcf" ;;
esac
