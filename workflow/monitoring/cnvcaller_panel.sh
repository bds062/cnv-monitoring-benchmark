#!/bin/bash
# CNVcaller discovery and genotyping on one panel, with the trio-benchmark settings
# (-f 0.1 -h 3 -r 0.5), from the moths' normalized read-depth files. CNVC_F overrides -f, and CNVC_OUT the
# monitoring folder written to (default $M): the -f 0 comparison sets CNVC_F=0 and CNVC_OUT to the monitoring
# folder of a second results tree ('f0' stage of run_all.sh; analysis/cnvcaller_f0.py).
#   bash cnvcaller_panel.sh PANEL                 -> $M/panels/PANEL/cnvcaller/PANEL.vcf.gz
#   bash cnvcaller_panel.sh PANEL SAMPLE [chr7]   the panel plus one new moth -> $M/query/PANEL/cnvcaller/SAMPLE.vcf.gz;
#                                                 chr7: read depth of chromosome 7 only (carrier titration)
source "$(dirname "$0")/common.sh"
P=$1; S=${2:-}; REGION=${3:-}
rd() { rd_of "$1"; }
O=${CNVC_OUT:-$M}
if [ -z "$S" ]; then W=$O/panels/$P/cnvcaller; NAME=$P; else W=$O/query/$P/cnvcaller/$S; NAME=$S; fi
mkdir -p "$W"; cd "$W"
cp "$CNVC_ZEA/referenceDB.800" .; : > exclude.txt; : > list.txt
for m in $(panel_members "$P") $S; do
    if [ "$REGION" == chr7 ]; then        # keep the file name: CNVcaller reads the sample name from it
        f=$(rd "$m"); mkdir -p rd; awk -v c="$CHR7" '$1 == c' "$f" > "rd/$(basename "$f")"; echo "$W/rd/$(basename "$f")"
    else rd "$m"; fi >> list.txt
done
sh "$CNVCALLER/CNV.Discovery.sh" -l list.txt -e exclude.txt -f "${CNVC_F:-0.1}" -h 3 -r 0.5 -p primaryCNVR -m mergeCNVR
python3 "$CNVCALLER/Genotype.py" --cnvfile mergeCNVR --outprefix "$NAME" --nproc "$THREADS"
bgzip -f "$NAME.vcf"
[ -z "$S" ] || mv "$NAME.vcf.gz" ..
