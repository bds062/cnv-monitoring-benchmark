#!/bin/bash
# align.sh SAMPLE [taylor2021]: Trimmomatic -> Bowtie2 --very-sensitive -> sorted, indexed BAM. Env: preprocess.
# Trimming regimes, as in the paper's runs:
#   trio        Trimmomatic 0.38, TruSeq3-PE.fa:2:30:10, SLIDINGWINDOW:5:20, MINLEN:50 (trio pipeline)
#   taylor2021  Trimmomatic 0.39, adapt_seq.fa:2:30:10:8:true, SLIDINGWINDOW:4:15, MINLEN:50 (Taylor et al. 2021)
#   taylor2024  Trimmomatic 0.39, TruSeq2-PE.fa:2:30:10, LEADING:3, TRAILING:3, SLIDINGWINDOW:4:30, MINLEN:36
#               (Taylor et al. 2024). The paper's runs used the founders' reads as Taylor et al. trimmed them,
#               which are not public; this re-trims the public reads. TruSeq2-PE.fa is the adapter file kept
#               with those trimmed reads; LEADING/TRAILING are Trimmomatic's documented defaults.
# Read-depth BAM (bam/zea/, all callers' Bowtie2 input), by group (published_moths.tsv):
#   trio moths, LA2002 -> trio;  LA2012, LA2017 -> taylor2021;  Benzon -> taylor2024
# With the argument taylor2021, an LA2002 or Benzon moth gets its second BAM, in $T21/bam/zea/ (see lib.sh).
# Trio moths keep their trimmed reads for speedseq.sh (BWA-MEM) and are also aligned to the GA-R assembly
# GCA_022343045.1 (bam_gca/zea/). Trimmomatic's summary goes to trim/zea/SAMPLE[.taylor2021].trimmomatic.log.
source "$(dirname "$0")/../lib.sh"
S=$1
TRIMMOMATIC_DIR=${TRIMMOMATIC_DIR:-$CNVB_DATA/tools}
TRIO=0; samples zea | grep -x "$S" > /dev/null && TRIO=1
case $TRIO:$(group_of "$S") in
1:)              DEFAULT=trio ;;
0:LA2002)        DEFAULT=trio ;;
0:LA2012|0:LA2017) DEFAULT=taylor2021 ;;
0:Benzon)        DEFAULT=taylor2024 ;;
*) echo "align.sh: $S is neither a trio moth nor in published_moths.tsv" >&2; exit 1 ;;
esac
REGIME=${2:-$DEFAULT}
if [ "$REGIME" = "$DEFAULT" ]; then OUT=$CNVB_RESULTS/bam/zea; TAG=
elif [ "$REGIME" = taylor2021 ] && [[ $(group_of "$S") =~ ^(LA2002|Benzon)$ ]]; then OUT=$T21/bam/zea; TAG=.taylor2021
else echo "align.sh: no $REGIME alignment of $S in the paper's runs" >&2; exit 1; fi
case $REGIME in
trio)       V=0.38; CLIP=$TRIMMOMATIC_DIR/Trimmomatic-0.38/adapters/TruSeq3-PE.fa:2:30:10; QUAL="SLIDINGWINDOW:5:20 MINLEN:50" ;;
taylor2021) V=0.39; CLIP=$(cd "$(dirname "$0")" && pwd)/adapt_seq.fa:2:30:10:8:true; QUAL="SLIDINGWINDOW:4:15 MINLEN:50" ;;
taylor2024) V=0.39; CLIP=$TRIMMOMATIC_DIR/Trimmomatic-0.39/adapters/TruSeq2-PE.fa:2:30:10
            QUAL="LEADING:3 TRAILING:3 SLIDINGWINDOW:4:30 MINLEN:36" ;;
esac
T=$CNVB_RESULTS/work/trimmed/zea; P=$T/$S.$REGIME          # trimmed reads: $P_1.fq.gz, $P_2.fq.gz
mkdir -p "$T" "$CNVB_RESULTS/trim/zea"

java -jar "$TRIMMOMATIC_DIR/Trimmomatic-$V/trimmomatic-$V.jar" PE -threads "$THREADS" \
    "$CNVB_DATA/reads/zea/${S}_1.fastq.gz" "$CNVB_DATA/reads/zea/${S}_2.fastq.gz" \
    "${P}_1.fq.gz" "${P}_1.unpaired.fq.gz" "${P}_2.fq.gz" "${P}_2.unpaired.fq.gz" \
    "ILLUMINACLIP:$CLIP" $QUAL 2> "$CNVB_RESULTS/trim/zea/$S$TAG.trimmomatic.log"
rm "${P}"_[12].unpaired.fq.gz

align() {                       # align INDEX OUTDIR
    mkdir -p "$2"
    bowtie2 --very-sensitive -p "$THREADS" -x "$1" -1 "${P}_1.fq.gz" -2 "${P}_2.fq.gz" \
        --rg-id "$S" --rg "SM:$S" | samtools sort -@ "$THREADS" -T "$P" -o "$2/$S.bam"
    samtools index "$2/$S.bam"
}
align "$(dirname "$REF_ZEA")/bowtie2/zea" "$OUT"
if [ $TRIO = 1 ]; then
    align "$CNVB_DATA/ref/zea_gca/bowtie2/gca" "$CNVB_RESULTS/bam_gca/zea"
else
    rm "${P}"_[12].fq.gz        # published moths: speedseq.sh takes their reads back out of the Bowtie2 BAM
fi
