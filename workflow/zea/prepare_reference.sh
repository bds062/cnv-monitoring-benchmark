#!/bin/bash
# H. zea reference files: indexes, the Z->X renamed copy for GATK-gCNV and CNVkit, the repeat BED
# (soft-masked runs of chromosomes 1-31), and the Bowtie2 index of the GA-R assembly. Env: preprocess.
source "$(dirname "$0")/../lib.sh"
D=$(dirname "$REF_ZEA")

samtools faidx "$REF_ZEA"
mkdir -p "$D/bowtie2"
bowtie2-build --threads "$THREADS" "$REF_ZEA" "$D/bowtie2/zea"
bwa index "$REF_ZEA"                            # for SpeedSeq (BWA-MEM) and the Delly2 mappability map

# Chromosomes 1-30 keep their number, 31 (Z) becomes X: CNVkit infers sex from a chromosome named X.
awk '/^>/ {id = substr($1, 2)
           if (match($0, / chromosome [0-9]+,/)) {n = substr($0, RSTART + 12, RLENGTH - 13); print id "\t" (n == 31 ? "X" : n)}
           else if ($0 ~ /mitochondrion/) print id "\tMT"
           else print id "\tunplaced_scaffold_" id}' "$REF_ZEA" > "$D/x_names.tsv"
awk 'NR == FNR {m[$1] = $2; next} /^>/ {print ">" m[substr($1, 2)]; next} {print}' \
    "$D/x_names.tsv" "$REF_ZEA" > "$D/zea_x.fna"
samtools faidx "$D/zea_x.fna"
samtools dict "$D/zea_x.fna" -o "$D/zea_x.dict"

mkdir -p "$CNVB_DATA/repeats"
# Repeat BED: NCBI WindowMasker soft-masking (lowercase), chromosomes 1-31, runs merged across lines
awk 'BEGIN {OFS = "\t"; s = -1}
     function flush() {if (s >= 0) print c, s, e; s = -1}
     /^>/ {flush(); c = substr($1, 2); keep = ($0 ~ / chromosome [0-9]+,/); p = 0; next}
     keep {line = $0
           while (match(line, /[a-z]+/)) {
               a = p + RSTART - 1; b = a + RLENGTH
               if (a == e && s >= 0) e = b; else {flush(); s = a; e = b}
               p += RSTART + RLENGTH - 1; line = substr(line, RSTART + RLENGTH)}
           p += length(line)}
     END {flush()}' "$REF_ZEA" > "$CNVB_DATA/repeats/zea_softmask.bed"

# GA-R assembly (GCA_022343045.1), the earlier reference, for the trio coverage/GC/MAPQ comparison
G=$CNVB_DATA/ref/zea_gca
mkdir -p "$G/bowtie2"
bowtie2-build --threads "$THREADS" "$G/GCA_022343045.1_ASM2234304v1_genomic.fna" "$G/bowtie2/gca"
