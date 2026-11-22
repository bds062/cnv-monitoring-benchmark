#!/bin/bash
# cram_to_bam.sh SAMPLE: 1000 Genomes GRCh38 CRAM -> BAM (all reads) in bam/human/. The released alignments
# are BWA-MEM, so the breakpoint callers read the same BAM (linked into bam_bwa/human/). Each caller restricts
# itself (GATK-gCNV: chr1-22 and X; CNVcaller and CNVkit: no ALT/HLA contigs). Env: preprocess.
source "$(dirname "$0")/../lib.sh"
S=$1
mkdir -p "$CNVB_RESULTS/bam/human" "$CNVB_RESULTS/bam_bwa/human"
samtools view -@ "$THREADS" -b -T "$REF_HUMAN" -o "$CNVB_RESULTS/bam/human/$S.bam" \
    "$CNVB_DATA/reads/human/$S.final.cram"
samtools index -@ "$THREADS" "$CNVB_RESULTS/bam/human/$S.bam"
ln -sf "../../bam/human/$S.bam" "$CNVB_RESULTS/bam_bwa/human/$S.bam"
ln -sf "../../bam/human/$S.bam.bai" "$CNVB_RESULTS/bam_bwa/human/$S.bam.bai"
