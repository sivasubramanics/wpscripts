#!/bin/bash

print_log() {
    echo "$(date "+%Y-%m-%d %H:%M:%S") $1"
}

if [ $# -ne 3 ]; then
    echo "Usage: $0 <input.fa> <out_prefix> <threads>"
    exit 1
fi

INPUT=$1
OUTPUT=$2
THREADS=$3
print_log "INPUT: ${INPUT} DB: lsatv11 OUTPUT: ${OUTPUT}_vs_lsatv11.psl LOG: ${OUTPUT}_vs_lsatv11.psl.log"
gmap -d lsatv11 \
    -D /lustre/BIF/nobackup/selva001/work/references/from_ncbi/Lactuca_sativa_v11/ncbi_dataset/data/GCF_002870075.4/gmap_build/ \
    -B 5 \
    -t ${THREADS} \
    -f psl \
    ${INPUT} > ${OUTPUT}_vs_lsatv11.psl 2> ${OUTPUT}_vs_lsatv11.psl.log
print_log "INPUT: ${INPUT} DB: lsat_pg OUTPUT: ${OUTPUT}_vs_lsat_pg.psl LOG: ${OUTPUT}_vs_lsat_pg.psl.log"
gmap -d lsat_pg \
    -D /lustre/BIF/nobackup/selva001/work/references/pangenomes/lsat/gmap_build/ \
    -B 5 \
    -t ${THREADS} \
    -f psl \
    ${INPUT} > ${OUTPUT}_vs_lsat_pg.psl 2> ${OUTPUT}_vs_lsat_pg.psl.log
print_log "INPUT: ${INPUT} DB: lsal_pg OUTPUT: ${OUTPUT}_vs_lsal_pg.psl LOG: ${OUTPUT}_vs_lsal_pg.psl.log"
gmap -d lsal_pg \
    -D /lustre/BIF/nobackup/selva001/work/references/pangenomes/lsal/gmap_build/ \
    -B 5 \
    -t ${THREADS} \
    -f psl \
    ${INPUT} > ${OUTPUT}_vs_lsal_pg.psl 2> ${OUTPUT}_vs_lsal_pg.psl.log

