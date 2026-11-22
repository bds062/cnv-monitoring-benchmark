# Settings shared by every workflow script (source it; do not run it).
# Paths default to folders inside the repository. To put data or results elsewhere,
# export CNVB_DATA / CNVB_RESULTS, or set them in config.local.sh (not tracked).
CNVB_REPO=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
[ -f "$CNVB_REPO/config.local.sh" ] && source "$CNVB_REPO/config.local.sh"

export CNVB_DATA=${CNVB_DATA:-$CNVB_REPO/data}
export CNVB_RESULTS=${CNVB_RESULTS:-$CNVB_REPO/results}
export THREADS=${THREADS:-8}

SAMPLES=$CNVB_REPO/data/samples
REF_ZEA=$CNVB_DATA/ref/zea/GCF_022581195.2_ilHelZeax1.1_genomic.fna
REF_HUMAN=$CNVB_DATA/ref/human/GRCh38_full_analysis_set_plus_decoy_hla.fa
CALLS=$CNVB_RESULTS/calls
