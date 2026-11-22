#!/bin/bash
# Download the public inputs into $CNVB_DATA (default: data/). See data/README.md.
#   bash data/download.sh [TARGET ...]    targets: zea_ref trios published human_ref crams truth repeats tools
# Everything runs in envs/preprocess.yml (sra-tools, pigz, ncbi-datasets-cli, curl, unzip); `tools` also needs
# git, make, a C++ compiler and apptainer.
set -euo pipefail
source "$(cd "$(dirname "$0")/.." && pwd)/config.sh"
DL=$CNVB_DATA

sra() {                         # sra RUN NAME: paired reads of one SRA run -> reads/zea/NAME_{1,2}.fastq.gz
    local d=$DL/reads/zea
    mkdir -p "$d/sra"
    prefetch --max-size 200G "$1" -O "$d/sra"
    fasterq-dump --split-files -e "$THREADS" -O "$d/sra" "$d/sra/$1"
    for r in 1 2; do pigz -p "$THREADS" -c "$d/sra/$1_$r.fastq" > "$d/${2}_$r.fastq.gz"; done
    rm -r "$d/sra/$1" "$d/sra/$1"_[12].fastq
}

get() { mkdir -p "$1"; curl -fsSL --retry 5 -C - -o "$1/$(basename "$2")" "$2"; }

target() {
    case $1 in
    zea_ref)        # NCBI RefSeq GCF_022581195.2 (ilHelZeax1.1), soft-masked by NCBI WindowMasker,
                    # and the earlier GA-R assembly GCA_022343045.1 (trio coverage/GC/MAPQ comparison)
        for a in zea:GCF_022581195.2 zea_gca:GCA_022343045.1; do
            mkdir -p "$DL/ref/${a%%:*}" && cd "$DL/ref/${a%%:*}"
            datasets download genome accession "${a#*:}" --include genome --filename "${a#*:}.zip"
            unzip -o -j "${a#*:}.zip" '*_genomic.fna' && rm "${a#*:}.zip"
        done ;;
    trios)          # parents: Taylor et al. 2024 (PRJNA1126687); children: this study
        awk -F'\t' 'NR > 1 {print $1, $3}' "$SAMPLES/trio_runs.tsv" | while read -r s run; do
            if [[ $run == SRR* ]]; then sra "$run" "$s"
            else echo "skip $s: reads released in the SRA with the paper, accession not yet public" >&2; fi
        done ;;
    published)      # 35 Louisiana moths (PRJNA751583) and 6 Benzon founders (PRJNA1126687)
        awk -F'\t' 'NR > 1 {print $1}' "$SAMPLES/published_moths.tsv" | while read -r run; do sra "$run" "$run"; done ;;
    human_ref)      # GRCh38 analysis set used by the 1000 Genomes high-coverage alignments
        get "$DL/ref/human" https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa ;;
    crams)          # 12 high-coverage 1000 Genomes CRAMs (Byrska-Bishop et al. 2022), md5-checked
        awk -F'\t' 'NR > 1' "$SAMPLES/human_crams.tsv" | while read -r s url md5; do
            mkdir -p "$DL/reads/human"
            curl -fsSL --retry 5 -C - -o "$DL/reads/human/$s.final.cram" "$url"
            curl -fsSL --retry 5 -o "$DL/reads/human/$s.final.cram.crai" "$url.crai"
            echo "$md5  $DL/reads/human/$s.final.cram" | md5sum -c -
        done ;;
    truth)          # long-read SV call sets for the human children
        get "$DL/truth/hgsvc2" https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HGSVC2/release/v2.0/integrated_callset/variants_freeze4_sv_insdel_sym.vcf.gz
        get "$DL/truth/hgsvc2" https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HGSVC2/release/v2.0/integrated_callset/variants_freeze4_sv_insdel.tsv.gz
        get "$DL/truth/hgsvc3" https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HGSVC3/release/Variant_Calls/1.0/GRCh38/variants_GRCh38_sv_insdel_sym_HGSVC2024v1.0.vcf.gz
        get "$DL/truth/platinum" https://platinum-pedigree-data.s3.amazonaws.com/truthset_v1.2/NA12878_hq_v1.2.svs.vcf.gz
        get "$DL/truth/platinum" https://platinum-pedigree-data.s3.amazonaws.com/truthset_v1.2/NA12878_hq_v1.2.svs.bed.gz ;;
    repeats)        # UCSC RepeatMasker track for hg38 (H. zea repeats come from the reference soft-masking)
        mkdir -p "$DL/repeats"
        curl -fsSL --retry 5 -o "$DL/repeats/hg38_rmsk.txt.gz" https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/rmsk.txt.gz ;;
    tools)          # caller software that conda does not provide
        local T=$DL/tools; mkdir -p "$T" && cd "$T"
        curl -fsSLO https://github.com/broadinstitute/gatk/releases/download/4.6.2.0/gatk-4.6.2.0.zip
        unzip -oq gatk-4.6.2.0.zip && rm gatk-4.6.2.0.zip
        for z in 5854850/Trimmomatic-0.38 5854859/Trimmomatic-0.39; do
            curl -fsSLO "https://github.com/usadellab/Trimmomatic/files/$z.zip"
            unzip -oq "${z#*/}.zip" && rm "${z#*/}.zip"
        done
        git clone https://github.com/JiangYuLab/CNVcaller.git
        git -C CNVcaller checkout d851794
        git -C CNVcaller apply "$CNVB_REPO/workflow/callers/cnvcaller.patch"
        git clone --recursive https://github.com/hall-lab/speedseq.git
        git -C speedseq checkout 4e60002 && git -C speedseq submodule update --init --recursive
        make -C speedseq
        apptainer pull smoove_0.2.8.sif docker://brentp/smoove:v0.2.8 ;;
    *) echo "unknown target $1" >&2; exit 1 ;;
    esac
}

for t in ${@:-zea_ref trios published human_ref crams truth repeats tools}; do target "$t"; done
