#!/bin/bash
# GATK-gCNV 4.6.2.0, 1 kb bins, --max-copy-number 10 in every model. Env: gatk.
#   gatk.sh prepare SPECIES                  bins, contig ploidy priors, interval shards
#   gatk.sh counts SPECIES SAMPLE            read counts
#   gatk.sh ploidy SPECIES BACKGROUND GROUP  contig ploidy, COHORT mode, on the background of GROUP
#   gatk.sh model SPECIES BACKGROUND GROUP SHARD   denoising/calling model, COHORT mode, one interval shard
#   gatk.sh call SPECIES BACKGROUND SAMPLE...  CASE-mode ploidy and calls against the background, postprocessed
# GROUP is the query trio, or "all" for backgrounds shared by every sample (backgrounds.tsv).
source "$(dirname "$0")/common.sh"
GATK=${GATK:-$CNVB_DATA/tools/gatk-4.6.2.0/gatk}
species_settings "$2"
G=$WORK/gatk/$SPECIES
export OMP_NUM_THREADS=$THREADS
counts() { for s in "$@"; do echo "-I $G/counts/$s.counts.hdf5"; done; }
shards() { printf '%02d\n' $(seq 0 $((GATK_SHARDS - 1))); }

case $1 in
prepare)
    mkdir -p "$G/shards"
    L=(); [ "$GATK_CONTIGS" = all ] || for c in $GATK_CONTIGS; do L+=(-L "$c"); done
    "$GATK" PreprocessIntervals -R "$REF_DEPTH" "${L[@]}" --bin-length 1000 \
        --interval-merging-rule OVERLAPPING_ONLY -O "$G/bins.interval_list"
    # priors: sex chromosome 1 or 2 copies, mitochondrion 1, everything else 2
    grep -v '^@' "$G/bins.interval_list" | cut -f1 | uniq \
        | awk -v OFS='\t' 'BEGIN {print "CONTIG_NAME", "PLOIDY_PRIOR_0", "PLOIDY_PRIOR_1", "PLOIDY_PRIOR_2", "PLOIDY_PRIOR_3"}
              $1 ~ /^(chr)?X$/ {print $1, 0.01, 0.49, 0.49, 0.01; next}
              $1 == "MT" {print $1, 0.01, 0.98, 0.01, 0; next}
              {print $1, 0.01, 0.01, 0.97, 0.01}' > "$G/ploidy_priors.tsv"
    # contiguous shards in genome order (whole-genome human models do not finish unsharded)
    grep -v '^@' "$G/bins.interval_list" > "$G/shards/body"
    split -d -a 2 -l $(( ($(wc -l < "$G/shards/body") + GATK_SHARDS - 1) / GATK_SHARDS )) "$G/shards/body" "$G/shards/part_"
    for i in $(shards); do
        cat <(grep '^@' "$G/bins.interval_list") "$G/shards/part_$i" > "$G/shards/shard_$i.interval_list"
    done
    rm "$G/shards/body" "$G/shards/part_"* ;;
counts)
    S=$3; mkdir -p "$G/counts"
    if [ "$SPECIES" = zea ]; then BAM=$(bam_x "$S" "$G/counts"); else BAM=$CNVB_RESULTS/bam/human/$S.bam; fi
    "$GATK" CollectReadCounts -I "$BAM" -L "$G/bins.interval_list" \
        --interval-merging-rule OVERLAPPING_ONLY -O "$G/counts/$S.counts.hdf5"
    [ "$SPECIES" = human ] || rm "$BAM" "$BAM.bai" ;;
ploidy)
    M=$G/$3/$4; mkdir -p "$M"
    export PYTENSOR_FLAGS=base_compiledir=$M/.pytensor_ploidy   # one compile cache per job: shared ones lock
    "$GATK" DetermineGermlineContigPloidy -L "$G/bins.interval_list" --interval-merging-rule OVERLAPPING_ONLY \
        --contig-ploidy-priors "$G/ploidy_priors.tsv" $(counts $(members "$SPECIES" "$3" "$4")) \
        --output "$M" --output-prefix ploidy ;;
model)
    M=$G/$3/$4; I=$5
    export PYTENSOR_FLAGS=base_compiledir=$M/.pytensor_$I
    "$GATK" GermlineCNVCaller --run-mode COHORT -L "$G/shards/shard_$I.interval_list" \
        --interval-merging-rule OVERLAPPING_ONLY --contig-ploidy-calls "$M/ploidy-calls" \
        $(counts $(members "$SPECIES" "$3" "$4")) --max-copy-number 10 \
        --output "$M" --output-prefix "shard_$I" ;;
call)
    # samples are CASE-called together (H. zea: one sample per call; human: the held-out trio)
    BG=$3; shift 3; M=$G/$BG/$(query_group "$SPECIES" "$BG" "$1"); W=$M/case_$(IFS=_; echo "$*")
    mkdir -p "$W"
    export PYTENSOR_FLAGS=base_compiledir=$W/.pytensor
    "$GATK" DetermineGermlineContigPloidy --model "$M/ploidy-model" $(counts "$@") \
        --output "$W" --output-prefix ploidy
    POST=()
    for i in $(shards); do
        "$GATK" GermlineCNVCaller --run-mode CASE --model "$M/shard_$i-model" \
            --contig-ploidy-calls "$W/ploidy-calls" $(counts "$@") --max-copy-number 10 \
            --output "$W" --output-prefix "shard_$i"
        POST+=(--calls-shard-path "$W/shard_$i-calls" --model-shard-path "$M/shard_$i-model")
    done
    for S in "$@"; do
        X=$(grep -lx "$S" "$W"/shard_00-calls/SAMPLE_*/sample_name.txt); X=${X%/sample_name.txt}
        O=$CALLS/$SPECIES/gatk/$BG/$S; mkdir -p "$O"
        "$GATK" --java-options -Xmx40g PostprocessGermlineCNVCalls "${POST[@]}" --sample-index "${X##*SAMPLE_}" \
            --contig-ploidy-calls "$W/ploidy-calls" \
            --output-genotyped-intervals "$O/$S.intervals.vcf.gz" \
            --output-genotyped-segments "$O/$S.segments.vcf.gz" \
            --output-denoised-copy-ratios "$O/$S.denoised_copy_ratios.tsv"
    done ;;
esac
