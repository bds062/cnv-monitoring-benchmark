# Species settings and sample-sheet lookups shared by the caller scripts (source it).
source "$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)/lib.sh"

WORK=$CNVB_RESULTS/work

species_settings() {            # species_settings SPECIES: reference and caller settings
    SPECIES=$1
    case $1 in
    zea)
        REF=$REF_ZEA
        REF_DEPTH=$CNVB_DATA/ref/zea/zea_x.fna      # GATK-gCNV and CNVkit: chromosome 31 (Z) named X
        X_NAMES=$CNVB_DATA/ref/zea/x_names.tsv      # GCF contig -> name in zea_x.fna
        SEX_CHR=NC_061482.1
        GATK_CONTIGS=all
        GATK_SHARDS=1
        PARX=() ;;
    human)
        REF=$REF_HUMAN
        REF_DEPTH=$CNVB_DATA/ref/human/GRCh38_noalt_decoy.fa   # CNVcaller, GATK-gCNV, CNVkit: no ALT/HLA contigs
        SEX_CHR=chrX
        GATK_CONTIGS="$(echo chr{1..22} chrX)"      # chrY left out of GATK-gCNV, as in the paper's runs
        GATK_SHARDS=30
        PARX=(--diploid-parx-genome grch38) ;;
    *) echo "unknown species $1" >&2; exit 1 ;;
    esac
}

trio_of() {                     # trio_of SPECIES SAMPLE
    awk -F'\t' -v s="$1" -v x="$2" '$1 == s && $4 == x {print $2}' "$SAMPLES/trios.tsv"
}

members() {                     # members SPECIES BACKGROUND TRIO: background of a query trio
    awk -F'\t' -v s="$1" -v b="$2" -v t="$3" '$1 == s && $2 == b && ($3 == t || $3 == "all") {print $4}' \
        "$SAMPLES/backgrounds.tsv"
}

query_group() {                 # query_group SPECIES BACKGROUND SAMPLE: TRIO, or "all" for shared backgrounds
    if awk -F'\t' -v s="$1" -v b="$2" '$1 == s && $2 == b && $3 == "all" {f = 1} END {exit !f}' "$SAMPLES/backgrounds.tsv"
    then echo all; else trio_of "$1" "$3"; fi
}

bam_x() {                       # bam_x SAMPLE DIR: H. zea Bowtie2 BAM with chromosome 31 named X (temporary copy)
    local out=$2/$1.x.bam
    samtools view -H "$CNVB_RESULTS/bam/zea/$1.bam" \
        | awk -F'\t' -v OFS='\t' 'NR == FNR {m[$1] = $2; next}
              /^@SQ/ {for (i = 2; i <= NF; i++) if ($i ~ /^SN:/) $i = "SN:" m[substr($i, 4)]} {print}' "$X_NAMES" - \
        > "$out.header.sam"
    samtools reheader "$out.header.sam" "$CNVB_RESULTS/bam/zea/$1.bam" > "$out"
    samtools index "$out"
    rm "$out.header.sam"
    echo "$out"
}

job() {                         # job ENV COMMAND...: run (or submit) COMMAND inside conda env ENV (envs/ENV.yml)
    local e=$1; shift
    run "conda run --no-capture-output -n $e $*"
}

groups() {                      # groups SPECIES BACKGROUND: query trios of a background ("all" if shared)
    awk -F'\t' -v s="$1" -v b="$2" '$1 == s && $2 == b {print $3}' "$SAMPLES/backgrounds.tsv" | sort -u
}
