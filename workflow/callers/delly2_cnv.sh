#!/bin/bash
# Delly2 read-depth copy number (delly cnv, 10 kb windows) on the BWA-MEM BAMs. Env: delly2.
#   delly2_cnv.sh prepare SPECIES        mappability map (dicey, 101 bp read pairs, as in the Delly docs)
#   delly2_cnv.sh call SPECIES SAMPLE    {sample}.vcf.gz (CNV calls) and {sample}.cov.gz (windows)
source "$(dirname "$0")/common.sh"
species_settings "$2"
M=$WORK/delly2_cnv/$SPECIES
case $1 in
prepare)
    mkdir -p "$M" && cd "$M"
    dicey chop "$REF"
    bwa mem -t "$THREADS" "$REF" read1.fq.gz read2.fq.gz | samtools sort -@ "$THREADS" -o srt.bam -
    samtools index srt.bam
    dicey mappability2 srt.bam
    gunzip -f map.fa.gz && bgzip map.fa && samtools faidx map.fa.gz
    rm read1.fq.gz read2.fq.gz srt.bam srt.bam.bai ;;
call)
    S=$3; O=$CALLS/$SPECIES/delly2_cnv/none/$S; mkdir -p "$O"
    delly cnv -a -g "$REF" -m "$M/map.fa.gz" -c "$O/$S.cov.gz" -o "$O/$S.bcf" "$CNVB_RESULTS/bam_bwa/$SPECIES/$S.bam"
    bcftools view "$O/$S.bcf" -Oz -o "$O/$S.vcf.gz"
    rm "$O/$S.bcf" "$O/$S.bcf.csi" ;;
esac
