# Settings for the chromosome 7 depth series (source it; do not run it).
source "$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)/monitoring/common.sh"

DS=$CNVB_RESULTS/downsample                   # everything the depth series writes
LEVELS=(full 10x 5x)                          # mean chromosome 7 depth of the Bowtie2 BAM; full = not subsampled
declare -A TARGET=([full]=0 [10x]=10 [5x]=5)
SEED=20261027
CALLERS=$(cd "$(dirname "${BASH_SOURCE[0]}")/../callers" && pwd)
# Each level is a results tree of its own (DS/LEVEL/{bam,bam_bwa,calls,work}), so the caller scripts of
# workflow/callers run on it unchanged with CNVB_RESULTS=DS/LEVEL (config.local.sh must not override it).
REF7=$M/chr7/ref/chr7.fa; MAP7=$M/chr7/ref/chr7.map.fa.gz
