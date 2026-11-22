#!/bin/bash
# GRCh38 reference files. The read-depth callers use the 1000 Genomes analysis set without its ALT and
# HLA contigs (primary assembly + decoys + EBV), so that windows with an ALT haplotype are not marked as
# duplicated. Repeat BED: UCSC RepeatMasker, chr1-22, X, Y, merged. Env: preprocess.
source "$(dirname "$0")/../lib.sh"
D=$(dirname "$REF_HUMAN")

samtools faidx "$REF_HUMAN"
samtools faidx "$REF_HUMAN" $(cut -f1 "$REF_HUMAN.fai" | grep -v -e '^HLA-' -e '_alt$') > "$D/GRCh38_noalt_decoy.fa"
samtools faidx "$D/GRCh38_noalt_decoy.fa"
samtools dict "$D/GRCh38_noalt_decoy.fa" -o "$D/GRCh38_noalt_decoy.dict"

zcat "$CNVB_DATA/repeats/hg38_rmsk.txt.gz" | awk -v OFS='\t' '$6 ~ /^chr([0-9]+|X|Y)$/ {print $6, $7, $8}' \
    | sort -k1,1 -k2,2n | bedtools merge > "$CNVB_DATA/repeats/human_rmsk.bed"
