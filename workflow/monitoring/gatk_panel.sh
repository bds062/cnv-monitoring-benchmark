#!/bin/bash
# GATK-gCNV on one panel, with the denoising and calling model restricted to chromosome 7.
#   bash gatk_panel.sh PANEL          COHORT mode: genome-wide ploidy, chromosome 7 model, every member
#                                     postprocessed -> $M/panels/PANEL/gatk/{member}.segments.vcf.gz
#   bash gatk_panel.sh PANEL SAMPLE   CASE mode against the panel's models (run the COHORT step first)
#                                     -> $M/query/PANEL/gatk/SAMPLE.segments.vcf.gz
source "$(dirname "$0")/common.sh"
P=$1; S=${2:-}; G=$M/panels/$P/gatk
counts() { while read -r m; do counts_of "$m"; done; }
# Each run gets its own cache and PyTensor compile directory: gCNV runs that share them (e.g. two CASE runs
# started in parallel inside one job) race on the arviz cache and the compile lock and fail. Run CASE moths
# one per job, or one after another within a job.
export XDG_CACHE_HOME=$(mktemp -d)
export PYTENSOR_FLAGS=base_compiledir=$XDG_CACHE_HOME; trap 'rm -rf "$XDG_CACHE_HOME"' EXIT

postprocess() {                # postprocess CALLS MODEL PLOIDY_CALLS INDEX OUT_PREFIX
    "$GATK" PostprocessGermlineCNVCalls --calls-shard-path "$1" --model-shard-path "$2" \
        --contig-ploidy-calls "$3" --sample-index "$4" \
        --output-genotyped-intervals "$5.intervals.vcf.gz" --output-genotyped-segments "$5.segments.vcf.gz" \
        --output-denoised-copy-ratios "$5.denoised.tsv"
}

if [ -z "$S" ]; then
    mkdir -p "$G"; cd "$G"
    species_settings zea                  # GATK-gCNV contig names: chromosome 7 is renamed in $X_NAMES
    C7=$(awk -v c="$CHR7" '$1 == c {print $2}' "$X_NAMES")
    awk -v c="$C7" '/^@/ || $1 == c' "$GATK_ZEA/bins.interval_list" > chr7.interval_list
    panel_members "$P" > members.txt
    mapfile -t IN < <(counts < members.txt | sed 's/^/--input\n/')
    "$GATK" DetermineGermlineContigPloidy -L "$GATK_ZEA/bins.interval_list" --interval-merging-rule OVERLAPPING_ONLY \
        --contig-ploidy-priors "$GATK_ZEA/ploidy_priors.tsv" "${IN[@]}" --output . --output-prefix ploidy
    "$GATK" GermlineCNVCaller --run-mode COHORT -L chr7.interval_list --interval-merging-rule OVERLAPPING_ONLY \
        --contig-ploidy-calls ploidy-calls "${IN[@]}" --max-copy-number 10 --output . --output-prefix bg
    i=0
    while read -r m; do postprocess bg-calls bg-model ploidy-calls $i "$m"; i=$((i + 1)); done < members.txt
else
    Q=$M/query/$P/gatk; mkdir -p "$Q/$S"; cd "$Q/$S"
    C=$(echo "$S" | counts)
    "$GATK" DetermineGermlineContigPloidy --model "$G/ploidy-model" --input "$C" \
        --interval-merging-rule OVERLAPPING_ONLY --output . --output-prefix ploidy
    "$GATK" GermlineCNVCaller --run-mode CASE --model "$G/bg-model" --contig-ploidy-calls ploidy-calls \
        --input "$C" --max-copy-number 10 --output . --output-prefix case
    postprocess case-calls "$G/bg-model" ploidy-calls 0 "$Q/$S"
fi
