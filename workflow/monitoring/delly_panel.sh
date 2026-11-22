#!/bin/bash
# Delly2 germline CNV workflow on one panel (Delly README, "Germline CNV calling"), on the
# chromosome 7 realignments of delly_chr7.sh: merge the moths' CNV calls into one site list
# (1 kb-1 Mb), genotype every moth at every site, merge the genotypes and filter them.
# MINQUAL is delly merge's minimum site quality: 0 (headline) or 200 (Delly2's default).
# QUERY moths are genotyped at the panel's sites without entering the merge or the filter,
# which is how a new moth is added to an existing panel (fixed query).
#   bash delly_panel.sh PANEL MINQUAL [QUERY...]
#   -> $M/panels/PANEL/delly_yMINQUAL/{sites,filtered}.bcf, $M/query/PANEL/delly_yMINQUAL/QUERY.bcf
source "$(dirname "$0")/common.sh"
P=$1; Y=$2; shift 2
REF=$M/chr7/ref/chr7.fa; MAP=$M/chr7/ref/chr7.map.fa.gz
W=$M/panels/$P/delly_y$Y; Q=$M/query/$P/delly_y$Y
mkdir -p "$W"; cd "$W"

genotype() {                    # genotype OUT.bcf SAMPLE
    delly cnv -a -u -v "$W/sites.bcf" -g "$REF" -m "$MAP" -o "$1" "$M/chr7/$2/$2.bam" && bcftools index -f "$1"
}

panel_members "$P" > members.txt
sed "s#.*#$M/chr7/&/&.cnvl.bcf#" members.txt > discovery.list
delly merge -e -p -m 1000 -n 1000000 -y "$Y" -o sites.bcf discovery.list && bcftools index -f sites.bcf
while read -r S; do genotype "geno_$S.bcf" "$S"; done < members.txt
sed "s#.*#geno_&.bcf#" members.txt > geno.list
bcftools merge -m id -O b -o merged.bcf -l geno.list && bcftools index -f merged.bcf
delly classify -f germline -o filtered.bcf merged.bcf && bcftools index -f filtered.bcf
for S in "$@"; do mkdir -p "$Q"; genotype "$Q/$S.bcf" "$S"; done
