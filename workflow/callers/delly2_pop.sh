#!/bin/bash
# Delly2 population read-depth workflow (the germline CNV workflow of the Delly README), the copy number of Delly2
# at Tryp77 and the Delly2 workflow of the monitoring analyses. Env: delly2.
#   delly2_pop.sh cnvl SPECIES SAMPLE              delly cnv -l: read-depth CNVs of one moth, segment ends
#                                                  refined by its own breakpoints (delly2.sh discovery BCF)
#   delly2_pop.sh panel SPECIES BACKGROUND SAMPLE|all
#       the members (background moths plus the query, as in delly2_merge.sh): delly merge -e -p -m 1000
#       -n 1000000 -y MERGE_Y of their delly cnv -l calls -> (chromosome 7 sites only if POP_CHR7=1) ->
#       delly cnv -u at the merged sites for every member (sample renamed to the member) -> bcftools merge ->
#       delly classify -f germline -> $CALLS/SPECIES/delly2_pop/BACKGROUND/SAMPLE/SAMPLE.bcf (every member's
#       column; with "all", one copy per member)
# Settings (environment): MERGE_Y (default 0, as in the monitoring analyses; Delly2's default 200 keeps no site
#   at Tryp77), POP_CHR7 (1: genotype the chromosome 7 sites only, enough for copy number at Tryp77 and
#   identical there to genotyping every site), BREAK_REF / MAP (depth series: chr7.fa and its map).
source "$(dirname "$0")/common.sh"
species_settings "$2"
R=${BREAK_REF:-$REF}; MP=${MAP:-$WORK/delly2_cnv/$SPECIES/map.fa.gz}
case $1 in
cnvl)
    S=$3; O=$CALLS/$SPECIES/delly2_pop/cnvl/$S; mkdir -p "$O"
    delly cnv -a -l "$CALLS/$SPECIES/delly2/none/$S/$S.bcf" -g "$R" -m "$MP" -o "$O/$S.bcf" \
        "$CNVB_RESULTS/bam_bwa/$SPECIES/$S.bam" && bcftools index -f "$O/$S.bcf" ;;
panel)
    BG=$3; S=$4; Y=${MERGE_Y:-0}
    W=$WORK/delly2_pop/$SPECIES/$BG/$S; mkdir -p "$W" && cd "$W"
    if [ "$S" = all ]; then members "$SPECIES" "$BG" all > members.txt
    else { members "$SPECIES" "$BG" "$(query_group "$SPECIES" "$BG" "$S")"; echo "$S"; } > members.txt; fi
    sed "s#.*#$CALLS/$SPECIES/delly2_pop/cnvl/&/&.bcf#" members.txt > discovery.list
    delly merge -e -p -m 1000 -n 1000000 -y "$Y" -o sites.bcf discovery.list && bcftools index -f sites.bcf
    if [ "${POP_CHR7:-0}" = 1 ]; then                                # H. zea chromosome 7
        bcftools view -r NC_061458.1 -Ob -o sites7.bcf sites.bcf && mv sites7.bcf sites.bcf
        bcftools index -f sites.bcf
    fi
    while read -r m; do
        delly cnv -a -u -v sites.bcf -g "$R" -m "$MP" -o "raw_$m.bcf" "$CNVB_RESULTS/bam_bwa/$SPECIES/$m.bam"
        echo "$m" > "name_$m.txt"                                   # every BAM names its sample L002
        bcftools reheader -s "name_$m.txt" -o "geno_$m.bcf" "raw_$m.bcf" && bcftools index -f "geno_$m.bcf"
    done < members.txt
    bcftools merge -m id -Ob -o merged.bcf $(sed 's#.*#geno_&.bcf#' members.txt) && bcftools index -f merged.bcf
    delly classify -f germline -o filtered.bcf merged.bcf && bcftools index -f filtered.bcf
    for m in $([ "$S" = all ] && cat members.txt || echo "$S"); do
        O=$CALLS/$SPECIES/delly2_pop/$BG/$m; mkdir -p "$O" && cp filtered.bcf "$O/$m.bcf" && cp filtered.bcf.csi "$O/$m.bcf.csi"
    done
    cd .. && rm -r "$W" ;;
esac
