#!/bin/bash
# delly2.sh SPECIES SAMPLE: Delly2 breakpoint discovery in one sample (delly call). The BCF is the input of
# delly2_merge.sh; the VCF is the sample's own calls (background "none"). BREAK_REF (environment) replaces the
# species reference (depth series: chr7.fa). Env: delly2.
source "$(dirname "$0")/common.sh"
species_settings "$1"; S=$2
O=$CALLS/$SPECIES/delly2/none/$S; mkdir -p "$O"
delly call -g "${BREAK_REF:-$REF}" -o "$O/$S.bcf" "$CNVB_RESULTS/bam_bwa/$SPECIES/$S.bam"
bcftools view "$O/$S.bcf" -Oz -o "$O/$S.vcf.gz"
