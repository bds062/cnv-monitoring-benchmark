#!/bin/bash
# Chromosome 7 depth series of the 12 trio moths: every caller at full depth, 10x and 5x, each with its
# headline background (crosstrio-9; CNVkit crosstrio-3) built from moths at the same depth.
#   chr7.sh prepare                  chromosome 7 parts of the GATK-gCNV bins, CNVkit 5 kb bins and CNVcaller
#                                    window database, staged into every level (once)
#   chr7.sh bams SAMPLE              chromosome 7 BAMs of one moth at every level -> DS/LEVEL/{bam,bam_bwa}/zea/
#   chr7.sh delly_cnv LEVEL SAMPLE   delly cnv on the BWA-MEM BAM (copy number at Tryp77)
#   chr7.sh gatk_counts LEVEL SAMPLE
#   chr7.sh gatk_model LEVEL TRIO    contig ploidy and COHORT model on the crosstrio-9 moths of TRIO
#   chr7.sh gatk_case LEVEL SAMPLE   CASE ploidy and calls against the model of the moth's trio, postprocessed
#   chr7.sh cnvkit_coverage LEVEL SAMPLE
#   chr7.sh cnvkit_call LEVEL SAMPLE reference from crosstrio-3, then fix, segment, call, export
# Delly2 and smoove discovery and merging and CNVcaller run through workflow/callers (see run_all.sh).
source "$(dirname "$0")/common.sh"
species_settings zea
C7=$(awk -v c="$CHR7" '$1 == c {print $2}' "$X_NAMES")         # chromosome 7 in the GATK-gCNV/CNVkit names
REFD=$DS/ref
L=${2:-}; S=${3:-}; T=$DS/$L

x_bam() {                       # x_bam IN OUT: chromosome 7 BAM with the GATK-gCNV/CNVkit contig name
    samtools view -H "$1" | sed "s/SN:$CHR7\b/SN:$C7/" > "$2.header.sam"
    samtools reheader "$2.header.sam" "$1" > "$2" && samtools index "$2" && rm -f "${2:?}.header.sam"
}
depth() { samtools coverage -r "$2" "$1" | awk 'NR == 2 {print $7}'; }     # depth BAM CONTIG: mean depth

case $1 in
prepare)
    mkdir -p "$REFD"
    awk -v c="$C7" '/^@/ || $1 == c' "$GATK_ZEA/bins.interval_list" > "$REFD/chr7.interval_list"
    printf 'CONTIG_NAME\tPLOIDY_PRIOR_0\tPLOIDY_PRIOR_1\tPLOIDY_PRIOR_2\tPLOIDY_PRIOR_3\n%s\t0.01\t0.01\t0.97\t0.01\n' "$C7" \
        > "$REFD/ploidy_priors.tsv"
    for b in target antitarget; do      # 5 kb bins of the matched protocol (WGS mode: the antitarget set is empty)
        awk -v c="$C7" '$1 == c' "$WORK/cnvkit/zea/bins5k/"*".$b.bed" > "$REFD/chr7.$b.bed"
    done
    # CNVcaller: windows and duplicated-window record of chromosome 7 only, so normalisation sees only chr7
    awk -v c="$CHR7" '$1 == c' "$CNVC_ZEA/referenceDB.800" > "$REFD/referenceDB.800"
    awk -v c="$CHR7" 'index($1, c ":") == 1' "$CNVC_ZEA/all_windows.link" > "$REFD/all_windows.link"
    for l in "${LEVELS[@]}"; do
        mkdir -p "$DS/$l/work/cnvcaller/zea" && cp "$REFD/referenceDB.800" "$REFD/all_windows.link" "$DS/$l/work/cnvcaller/zea/"
    done ;;
bams)
    S=$2; W=$(mktemp -d); trap 'rm -rf "${W:?}"' EXIT
    samtools view -b "$CNVB_RESULTS/bam/zea/$S.bam" "$CHR7" \
        | samtools addreplacerg -r "ID:$S" -r "SM:$S" -r "PL:ILLUMINA" -o "$W/bt2.bam" - && samtools index "$W/bt2.bam"
    BWA=$M/chr7/$S/$S.bam                # chromosome 7 BWA-MEM realignment of the monitoring workflow
    D0=$(depth "$W/bt2.bam" "$CHR7")
    mkdir -p "$DS/depth"; printf 'level\tfraction\tdepth_bt2\tdepth_bwa\tfull_depth_bt2\n' > "$DS/depth/$S.tsv"
    for l in "${LEVELS[@]}"; do
        mkdir -p "$DS/$l/bam/zea" "$DS/$l/bam_bwa/zea"; O=$DS/$l/bam/zea/$S.bam; OB=$DS/$l/bam_bwa/zea/$S.bam
        # fraction = target / the moth's Bowtie2 depth, and the same fraction and seed for its BWA-MEM BAM;
        # samtools -s keeps or drops whole read pairs. A moth already at or below the target is not subsampled.
        F=$(awk -v d="$D0" -v x="${TARGET[$l]}" 'BEGIN {printf "%.6f", (x == 0 || d <= x) ? 1 : x / d}')
        if [ "$F" = 1.000000 ]; then cp "$W/bt2.bam" "$O"; samtools view -b -o "$OB" "$BWA"
        else samtools view -b -s "$SEED${F#0}" -o "$O" "$W/bt2.bam"; samtools view -b -s "$SEED${F#0}" -o "$OB" "$BWA"; fi
        samtools index "$O" && samtools index "$OB"
        printf '%s\t%s\t%s\t%s\t%s\n' "$l" "$F" "$(depth "$O" "$CHR7")" "$(depth "$OB" "$CHR7")" "$D0" >> "$DS/depth/$S.tsv"
    done ;;
delly_cnv)
    O=$T/calls/zea/delly2_cnv/none/$S; mkdir -p "$O"
    delly cnv -a -g "$REF7" -m "$MAP7" -c "$O/$S.cov.gz" -o "$O/$S.bcf" "$T/bam_bwa/zea/$S.bam" ;;
gatk_counts)
    G=$T/work/gatk; mkdir -p "$G/counts"; x_bam "$T/bam/zea/$S.bam" "$G/counts/$S.x.bam"
    "$GATK" CollectReadCounts -I "$G/counts/$S.x.bam" -L "$REFD/chr7.interval_list" \
        --interval-merging-rule OVERLAPPING_ONLY -O "$G/counts/$S.counts.hdf5"
    rm -f "${G:?}/counts/${S:?}.x.bam" "${G:?}/counts/${S:?}.x.bam.bai" ;;
gatk_model)
    # contig ploidy is fitted on the chromosome 7 bins: it sets each moth's depth, which downsampling changes
    Q=$3; G=$T/work/gatk; W=$G/$Q; mkdir -p "$W"
    export PYTENSOR_FLAGS=base_compiledir=$W/.pytensor
    IN=(); for m in $(members zea crosstrio9 "$Q"); do IN+=(--input "$G/counts/$m.counts.hdf5"); done
    "$GATK" DetermineGermlineContigPloidy -L "$REFD/chr7.interval_list" --interval-merging-rule OVERLAPPING_ONLY \
        --contig-ploidy-priors "$REFD/ploidy_priors.tsv" "${IN[@]}" --output "$W" --output-prefix ploidy
    "$GATK" GermlineCNVCaller --run-mode COHORT -L "$REFD/chr7.interval_list" --interval-merging-rule OVERLAPPING_ONLY \
        --contig-ploidy-calls "$W/ploidy-calls" "${IN[@]}" --max-copy-number 10 --output "$W" --output-prefix bg ;;
gatk_case)
    G=$T/work/gatk; Q=$(trio_of zea "$S"); W=$G/case/$S; O=$T/calls/zea/gatk/crosstrio9/$S; mkdir -p "$W" "$O"
    export PYTENSOR_FLAGS=base_compiledir=$W/.pytensor
    "$GATK" DetermineGermlineContigPloidy --model "$G/$Q/ploidy-model" --input "$G/counts/$S.counts.hdf5" \
        --interval-merging-rule OVERLAPPING_ONLY --output "$W" --output-prefix ploidy
    "$GATK" GermlineCNVCaller --run-mode CASE --model "$G/$Q/bg-model" --contig-ploidy-calls "$W/ploidy-calls" \
        --input "$G/counts/$S.counts.hdf5" --max-copy-number 10 --output "$W" --output-prefix case
    "$GATK" PostprocessGermlineCNVCalls --calls-shard-path "$W/case-calls" --model-shard-path "$G/$Q/bg-model" \
        --contig-ploidy-calls "$W/ploidy-calls" --output-genotyped-intervals "$O/$S.intervals.vcf.gz" \
        --output-genotyped-segments "$O/$S.segments.vcf.gz" --output-denoised-copy-ratios "$O/$S.denoised_copy_ratios.tsv" ;;
cnvkit_coverage)
    K=$T/work/cnvkit; mkdir -p "$K"; x_bam "$T/bam/zea/$S.bam" "$K/$S.x.bam"
    cnvkit.py coverage "$K/$S.x.bam" "$REFD/chr7.target.bed" -p "$THREADS" -o "$K/$S.targetcoverage.cnn"
    cnvkit.py coverage "$K/$S.x.bam" "$REFD/chr7.antitarget.bed" -p 1 -o "$K/$S.antitargetcoverage.cnn"
    rm -f "${K:?}/${S:?}.x.bam" "${K:?}/${S:?}.x.bam.bai" ;;
cnvkit_call)
    K=$T/work/cnvkit; O=$T/calls/zea/cnvkit/crosstrio3/$S; mkdir -p "$O" && cd "$O"
    REFS=(); for m in $(members zea crosstrio3 "$(trio_of zea "$S")"); do
        REFS+=("$K/$m.targetcoverage.cnn" "$K/$m.antitargetcoverage.cnn"); done
    cnvkit.py reference "${REFS[@]}" -f "$REF_DEPTH" -o reference.cnn
    cnvkit.py fix "$K/$S.targetcoverage.cnn" "$K/$S.antitargetcoverage.cnn" reference.cnn -o "$S.cnr"
    cnvkit.py segment "$S.cnr" -p 4 -o "$S.cns"
    cnvkit.py call "$S.cns" -o "$S.call.cns"
    cnvkit.py export vcf "$S.call.cns" -i "$S" -o "$S.call.vcf" && bgzip -f "$S.call.vcf" ;;
*)
    sed -n '2,13p' "$0"; exit 1 ;;
esac
