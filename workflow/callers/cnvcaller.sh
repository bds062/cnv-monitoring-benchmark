#!/bin/bash
# CNVcaller (d851794 + cnvcaller.patch), 800 bp windows. Env: cnvcaller.
#   cnvcaller.sh prepare SPECIES               window database and duplicated-window record
#   cnvcaller.sh rd SPECIES SAMPLE             read depth of one sample (sex chromosome scaled separately)
#   cnvcaller.sh call SPECIES BACKGROUND SAMPLE  CNV regions and genotypes of the background plus this one sample
source "$(dirname "$0")/common.sh"
CNVCALLER=${CNVCALLER_DIR:-$CNVB_DATA/tools/CNVcaller}
species_settings "$2"
C=$WORK/cnvcaller/$SPECIES
# H. zea reads are on the GCF contig names; the human record is built on the analysis set without ALT/HLA
[ "$SPECIES" = zea ] && FA=$REF || FA=$REF_DEPTH

dup_windows() {                 # dup_windows FASTA OUT: 800 bp k-mers aligned back with blasr, multi-copy windows kept
    python3 "$CNVCALLER/bin/0.1.Kmer_Generate.py" "$1" 800 "$2.kmer.fa"
    sawriter "$2.sa" "$1"
    blasr "$2.kmer.fa" "$1" --sa "$2.sa" --out "$2.aln" -m 5 --noSplitSubreads --minMatch 15 --maxMatch 20 \
        --advanceHalf --advanceExactMatches 10 --fastMaxInterval --fastSDP --aggressiveIntervalCut --bestn 10 \
        --nproc "$THREADS"
    python3 "$CNVCALLER/bin/0.2.Kmer_Link.py" "$2.aln" 800 "$2"
    rm "$2.kmer.fa" "$2.aln" "$2.sa"
}

case $1 in
prepare)
    mkdir -p "$C" && cd "$C"
    perl "$CNVCALLER/bin/CNVReferenceDB.pl" "$FA" -w 800
    if [ "$SPECIES" = zea ]; then       # H. zea: each sequence against itself, as in the paper's run
        mkdir -p split; : > all_windows.link
        for c in $(cut -f1 "$FA.fai"); do
            samtools faidx "$FA" "$c" > "split/$c.fa"
            dup_windows "split/$c.fa" "split/$c.link"
            cat "split/$c.link" >> all_windows.link
        done
    else                                # human: whole genome at once
        dup_windows "$FA" all_windows.link
    fi ;;
rd)
    S=$3; D=$C/rd/$S; rm -rf "$D"; mkdir -p "$D" && cd "$D"
    cp "$C/referenceDB.800" .
    sh "$CNVCALLER/Individual.Process.sh" -b "$CNVB_RESULTS/bam/$SPECIES/$S.bam" -h "$S" \
        -d "$C/all_windows.link" -s "$SEX_CHR" ;;
call)
    BG=$3; S=$4; W=$C/$BG/$S; O=$CALLS/$SPECIES/cnvcaller/$BG/$S
    rm -rf "$W"; mkdir -p "$W" "$O" && cd "$W"
    cp "$C/referenceDB.800" .
    for m in $(members "$SPECIES" "$BG" "$(query_group "$SPECIES" "$BG" "$S")") $S; do
        ls "$C/rd/$m/RD_normalized/"*
    done | awk '!seen[$0]++' > rd_list.txt
    : > exclude.txt
    sh "$CNVCALLER/CNV.Discovery.sh" -l rd_list.txt -e exclude.txt -f 0.1 -h 3 -r 0.5 -p primaryCNVR -m mergeCNVR
    python3 "$CNVCALLER/Genotype.py" --cnvfile mergeCNVR --outprefix "$S" --nproc 4
    bgzip -c "$S.vcf" > "$O/$S.vcf.gz" ;;
esac
