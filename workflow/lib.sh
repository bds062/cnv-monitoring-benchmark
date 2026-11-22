# Helpers for workflow scripts: source config.sh, fail early, and run one job either
# locally or through a scheduler. Set SUBMIT to submit instead of running, e.g.
#   SUBMIT="sbatch -c 8 --mem 32G -t 24:00:00 --wrap" bash workflow/zea/run_all.sh
set -euo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)/config.sh"

run() {                         # run "command string" now, or hand it to $SUBMIT
    if [ -n "${SUBMIT:-}" ]; then $SUBMIT "$*"; else bash -c "$*"; fi
}

samples() {                     # samples SPECIES [ROLE]: sample names from data/samples/trios.tsv
    awk -F'\t' -v s="$1" -v r="${2:-}" 'NR > 1 && $1 == s && (r == "" || $3 == r) {print $4}' "$SAMPLES/trios.tsv"
}

group_of() {                    # group_of SAMPLE: LA2002, LA2012, LA2017 or Benzon (published_moths.tsv); empty otherwise
    awk -F'\t' -v s="$1" 'NR > 1 && $1 == s {print $4}' "$SAMPLES/published_moths.tsv"
}

# The paper's runs aligned two groups of published moths twice (workflow/zea/align.sh): their read-depth BAM
# (bam/zea/) is trimmed like the trio moths (LA2002) or as Taylor et al. (2024) trimmed them (Benzon), and a
# second Bowtie2 BAM from Taylor et al. (2021)-style trimming sits in the tree $CNVB_RESULTS/t21. That second
# BAM feeds SpeedSeq (both groups) and the monitoring analyses (LA2002). Other moths have one BAM.
T21=$CNVB_RESULTS/t21
bam_taylor2021() {              # bam_taylor2021 SAMPLE: Bowtie2 BAM of a published moth from Taylor 2021-style trimming
    case $(group_of "$1") in
    LA2002|Benzon) echo "$T21/bam/zea/$1.bam" ;;
    *)             echo "$CNVB_RESULTS/bam/zea/$1.bam" ;;
    esac
}
