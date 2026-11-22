#!/bin/bash
# Caller-free truth for the published moths, from their Bowtie2 BAMs (Taylor 2021 trimming, common.sh):
#   depth     mean genome-wide depth (mapped reads x 150 bp / genome length)
#   depth_cn  read-depth copy number at the ddPCR amplicon: 2 x mean depth over the amplicon +-2 kb
#             (five 1 kb windows) / median depth over chromosome 7 (50 kb windows); carrier if >= 2.5
#   bash read_depth.sh   -> $M/read_depth.tsv
source "$(dirname "$0")/common.sh"
AMP_MID=8252641                               # ddPCR amplicon 8,252,593-8,252,690
mkdir -p "$M"; T=$(mktemp -d)
printf 'sample\tdepth\tdepth_cn\n' > "$M/read_depth.tsv"
for S in $(panel_members '2002+2012+2017'); do
    BAM=$(bt2_bam "$S")
    len=$(samtools idxstats "$BAM" | awk -v c="$CHR7" '$1 == c {print $2}')
    awk -v n="$len" -v c="$CHR7" 'BEGIN {for (a = 0; a < n - 50000; a += 50000) print c "\t" a "\t" a + 50000}' > "$T/chr7.bed"
    awk -v m="$AMP_MID" -v c="$CHR7" 'BEGIN {b = int(m / 1000) * 1000; for (a = b - 2000; a <= b + 2000; a += 1000) print c "\t" a "\t" a + 1000}' > "$T/amp.bed"
    med=$(samtools bedcov "$T/chr7.bed" "$BAM" | awk '{print $4 / ($3 - $2)}' | sort -g \
          | awk '{v[NR] = $1} END {print (NR % 2) ? v[(NR + 1) / 2] : (v[NR / 2] + v[NR / 2 + 1]) / 2}')
    amp=$(samtools bedcov "$T/amp.bed" "$BAM" | awk '{s += $4 / ($3 - $2)} END {print s / NR}')
    samtools idxstats "$BAM" | awk -v s="$S" -v a="$amp" -v m="$med" \
        '{L += $2; r += $3} END {printf "%s\t%.4f\t%.4f\n", s, r * 150 / L, 2 * a / m}' >> "$M/read_depth.tsv"
done
rm -r "$T"
