#!/bin/bash
# speedseq.sh SAMPLE: genome-wide BWA-MEM alignment with SpeedSeq align (BWA-MEM, SAMBLASTER duplicate marking,
# sambamba sort), plus the split-read and discordant-pair BAMs, in bam_bwa/zea/. These BAMs feed the breakpoint
# callers (trio moths and the published background moths) and delly cnv (trio moths).
# Input reads, as in the paper's runs:
#   trio moths       their trimmed reads (align.sh, trio regime)
#   published moths  the read pairs taken back out of their Taylor 2021-trimmed Bowtie2 BAM (lib.sh
#                    bam_taylor2021): primary records only, mate suffixes .1/.2 stripped so that collate pairs
#                    them (needed for BAMs made from fastq-dump --readids reads; harmless otherwise)
# Env: speedseq (SpeedSeq bin/ on PATH).
source "$(dirname "$0")/../lib.sh"
S=$1
OUT=$CNVB_RESULTS/bam_bwa/zea; mkdir -p "$OUT"
W=$CNVB_RESULTS/work/speedseq/$S; rm -rf "$W"; mkdir -p "$W"
if samples zea | grep -x "$S" > /dev/null; then
    R1=$CNVB_RESULTS/work/trimmed/zea/$S.trio_1.fq.gz; R2=${R1%_1.fq.gz}_2.fq.gz
else
    R1=$W/r1.fq.gz; R2=$W/r2.fq.gz
    samtools view -h "$(bam_taylor2021 "$S")" \
        | awk 'BEGIN {OFS = "\t"} !/^@/ && split($1, a, ".") == 3 {$1 = a[1] "." a[2]} {print}' \
        | samtools collate -@ 4 -u -O - "$W/coll" \
        | samtools fastq -@ 4 -n -F 0x900 -1 "$R1" -2 "$R2" -0 /dev/null -s /dev/null -
fi
speedseq align -t "$THREADS" -T "$W/sq" -R "@RG\tID:$S\tSM:$S\tLB:$S\tPL:ILLUMINA" -o "$OUT/$S" \
    "$REF_ZEA" "$R1" "$R2"
rm -r "$W"
samples zea | grep -x "$S" > /dev/null && rm "$R1" "$R2"
true
