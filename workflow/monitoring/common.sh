# Settings and helpers for the monitoring workflow (source it; do not run it).
source "$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)/callers/common.sh"

M=$CNVB_RESULTS/monitoring                    # everything the monitoring workflow writes
CHR7=NC_061458.1                              # chromosome 7 of GCF_022581195.2, holds Tryp77
PANELS=$SAMPLES/timeseries_panels.tsv
GATK=${GATK:-$CNVB_DATA/tools/gatk-4.6.2.0/gatk}
CNVCALLER=${CNVCALLER_DIR:-$CNVB_DATA/tools/CNVcaller}
# per-sample inputs made by workflow/callers for all 47 moths: GATK-gCNV bins, ploidy priors and read
# counts, CNVcaller read depth, and the Delly2 mappability map
GATK_ZEA=$WORK/gatk/zea
CNVC_ZEA=$WORK/cnvcaller/zea
DELLY_MAP=$WORK/delly2_cnv/zea/map.fa.gz
# The monitoring analyses read every Louisiana moth as Taylor et al. (2021) trimmed it. For LA2002 that is not
# the read-depth BAM of the trio benchmark but the second BAM in $T21 (workflow/lib.sh), so their read counts
# and read depth are made again there ('prepare' stage of run_all.sh).
bt2_bam() { bam_taylor2021 "$1"; }                       # bt2_bam SAMPLE: Bowtie2 BAM of the monitoring analyses
tree_of() { [ "$(group_of "$1")" = LA2002 ] && echo "$T21" || echo "$CNVB_RESULTS"; }
counts_of() { echo "$(tree_of "$1")/work/gatk/zea/counts/$1.counts.hdf5"; }   # counts_of SAMPLE: GATK-gCNV counts
rd_of() { ls "$(tree_of "$1")/work/cnvcaller/zea/rd/$1/RD_normalized/"*; }      # rd_of SAMPLE: CNVcaller read depth

panel_members() {              # panel_members PANEL; crosstrio9_TrioN = the nine moths of the other three trios
    case $1 in
        crosstrio9_*) awk -F'\t' -v t="${1#crosstrio9_}" 'NR > 1 && $1 == "zea" && $2 != t {print $4}' "$SAMPLES/trios.tsv" ;;
        *) awk -F'\t' -v p="$1" 'NR > 1 && $1 == p {print $2}' "$PANELS" ;;
    esac
}

panels() {                     # panels REGEX: panel names matching an extended regex, in sheet order
    awk -F'\t' -v r="^($1)\$" 'NR > 1 && $1 ~ r && !s[$1]++ {print $1}' "$PANELS"
}

SERIES='2002|2012|2017|2002[+]2012|2002[+]2017|2012[+]2017|2002[+]2012[+]2017|mixed13_.'
GROWING='2002|2002[+]2012|2002[+]2012[+]2017'
