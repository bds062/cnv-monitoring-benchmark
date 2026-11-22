#!/bin/bash
# delly2_merge.sh SPECIES BACKGROUND SAMPLE|all: Delly2 germline workflow (Delly README) for one trio member
# against one background, so Delly2 sees the same background as the read-depth callers. The members are the
# background moths plus the trio member, never a parent with its own child; with "all" (pooled12, a background
# of the trio moths themselves) the members are the background and every one of them is written.
#   delly merge -y MERGE_Y of the members' discovery calls (delly2.sh) -> delly call -v at the merged sites for
#   every member (sample renamed to the member) -> bcftools merge -> delly filter -f germline -> the trio
#   member's column -> $CALLS/SPECIES/delly2/BACKGROUND[SUFFIX]/SAMPLE/SAMPLE.vcf.gz
# The merged sites are kept: $CNVB_RESULTS/delly2_sites/SPECIES/BACKGROUND[SUFFIX]/SAMPLE/sites.bcf.
# Settings (environment):
#   MERGE_Y       delly merge minimum site quality; default 0 (no minimum), the setting of every Delly2 run in
#                 the paper, because Delly2's default of 200 removes the Tryp77 amplification from its
#                 read-depth workflow (delly2_pop.sh); the trio MIR is the same with either setting
#   MERGE_SUFFIX  appended to the background's directory; default "" for -y 0, "_y<MERGE_Y>" otherwise
#                 (MERGE_Y=200: the comparison with Delly2's default in BACKGROUND_y200).
#   BREAK_REF     reference FASTA of the BWA-MEM BAMs; default the species reference (depth series: chr7.fa)
# Env: delly2.
source "$(dirname "$0")/common.sh"
species_settings "$1"; BG=$2; S=$3
Y=${MERGE_Y:-0}
if [ "$Y" = 0 ]; then BGO=$BG${MERGE_SUFFIX-}; else BGO=$BG${MERGE_SUFFIX-_y$Y}; fi
R=${BREAK_REF:-$REF}
W=$WORK/delly2_merge/$SPECIES/$BGO/$S; mkdir -p "$W" && cd "$W"
if [ "$S" = all ]; then members "$SPECIES" "$BG" all > members.txt
else { members "$SPECIES" "$BG" "$(query_group "$SPECIES" "$BG" "$S")"; echo "$S"; } > members.txt; fi
sed "s#.*#$CALLS/$SPECIES/delly2/none/&/&.bcf#" members.txt > discovery.list
delly merge -y "$Y" -o sites.bcf discovery.list && bcftools index -f sites.bcf
while read -r m; do
    delly call -g "$R" -v sites.bcf -o "raw_$m.bcf" "$CNVB_RESULTS/bam_bwa/$SPECIES/$m.bam"
    echo "$m" > "name_$m.txt"
    bcftools reheader -s "name_$m.txt" -o "geno_$m.bcf" "raw_$m.bcf" && bcftools index -f "geno_$m.bcf"
done < members.txt
bcftools merge -m id -Ob -o merged.bcf $(sed 's#.*#geno_&.bcf#' members.txt) && bcftools index -f merged.bcf
delly filter -f germline -o filtered.bcf merged.bcf
for m in $([ "$S" = all ] && cat members.txt || echo "$S"); do
    O=$CALLS/$SPECIES/delly2/$BGO/$m; mkdir -p "$O"
    bcftools view -s "$m" filtered.bcf -Oz -o "$O/$m.vcf.gz"
done
K=$CNVB_RESULTS/delly2_sites/$SPECIES/$BGO/$S; mkdir -p "$K" && cp sites.bcf sites.bcf.csi "$K/"
cd .. && rm -r "$W"
