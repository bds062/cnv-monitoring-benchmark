#!/bin/bash
# Delly2 per moth on chromosome 7: realign the read pairs placed on chromosome 7 with BWA-MEM to
# chromosome 7 alone, remove duplicates, then delly call (breakpoints), delly cnv (read-depth copy
# number) and delly cnv -l (CNV calls refined with the moth's own breakpoints, the input of the
# germline merge workflow). Trio moths start from their BWA-MEM (SpeedSeq) BAMs, published moths
# from their Bowtie2 BAMs (Taylor 2021 trimming, common.sh).
#   bash delly_chr7.sh ref            chromosome 7 reference, BWA index and mappability slice (once)
#   bash delly_chr7.sh SAMPLE         one moth -> $M/chr7/SAMPLE/SAMPLE.{bam,bcf,cov.gz,cnvl.bcf}
source "$(dirname "$0")/common.sh"
REF=$M/chr7/ref/chr7.fa; MAP=$M/chr7/ref/chr7.map.fa.gz
BWA=${SPEEDSEQ_BWA:-$CNVB_DATA/tools/speedseq/bin/bwa}     # SpeedSeq's BWA 0.7.10, as in the paper's run

if [ "$1" == ref ]; then
    mkdir -p "$M/chr7/ref"
    samtools faidx "$REF_ZEA" "$CHR7" > "$REF" && samtools faidx "$REF" && "$BWA" index "$REF"
    samtools faidx "$DELLY_MAP" "$CHR7" | bgzip > "$MAP" && samtools faidx "$MAP"
    exit
fi

S=$1; W=$M/chr7/$S
if cut -f4 "$SAMPLES/trios.tsv" | grep -qx "$S"; then BAM=$CNVB_RESULTS/bam_bwa/zea/$S.bam; else BAM=$(bt2_bam "$S"); fi
mkdir -p "$W"; cd "$W"
samtools view -u -F 0x900 "$BAM" "$CHR7" | samtools collate -O -u - tmp \
    | samtools fastq -1 r1.fq.gz -2 r2.fq.gz -s /dev/null -0 /dev/null -
"$BWA" mem -t "$THREADS" -R "@RG\tID:$S\tSM:$S\tPL:ILLUMINA" "$REF" r1.fq.gz r2.fq.gz \
    | samtools fixmate -m -u - - | samtools sort -u - | samtools markdup -r - "$S.bam"
samtools index "$S.bam" && rm r1.fq.gz r2.fq.gz
delly call -g "$REF" -o "$S.bcf" "$S.bam"
delly cnv -a -g "$REF" -m "$MAP" -c "$S.cov.gz" -o "$S.cnv.bcf" "$S.bam"
delly cnv -a -l "$S.bcf" -g "$REF" -m "$MAP" -o "$S.cnvl.bcf" "$S.bam" && bcftools index -f "$S.cnvl.bcf"
