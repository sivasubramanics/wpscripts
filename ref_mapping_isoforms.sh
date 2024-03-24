#!/bin/bash
# script: ref_mapping_isoforms.sh
# description: map the isoform/transcripts on the genome using minimap2.
# contact: siva (siva.subramani975[at]gmail.com)
# date: 06-12-2023

set -e
set -u


# function to print log
function print_log(){
    echo "[`date +'%Y-%m-%d %H:%M:%S'`] LOG: $1"
}


# function to run command
function run_cmd(){
    echo "[`date +'%Y-%m-%d %H:%M:%S'`] CMD: $1"
    eval $1
}

# function to print usage
function usage(){
    echo "Usage: $0 [-h] [-f genome.fa] [-g genome.gtf] [-t transcripts.fa] [-o output_prefix] [-p num_threads]"
    echo "Options:"
    echo "  -f STR    genome fasta file"
    echo "  -g STR    genome gtf file"
    echo "  -t STR    transcripts fasta file"
    echo "  -o STR    output prefix"
    echo "  -p INT    number of threads"
    echo "  -h        help"
    exit 1
}

# function to convert seconds to hh:mm:ss
function seconds_to_hhmmss(){
    local T=$1
    local H=$((T/3600))
    local M=$((T%3600/60))
    local S=$((T%60))
    printf "%02d:%02d:%02d\n" $H $M $S
}

# check tools if installed
for tool in minimap2 bedToGenePred paftools.js gffcompare; do
    if ! command -v $tool &> /dev/null; then
        case $tool in
            minimap2)
                echo "minimap2 is not installed. You may install using the below command"
                echo "mamba install minimap2"
                exit 1
                ;;
            bedToGenePred)
                echo "bedToGenePred is not installed. You may install using the below command"
                echo "mamba install ucsc-bedtogenepred"
                exit 1
                ;;
            paftools.js)
                echo "paftools.js is not installed. You may install using the below command"
                echo "mamba install minimap2"
                exit 1
                ;;
            gffcompare)
                echo "gffcompare is not installed. You may install using the below command"
                echo "mamba install gffcompare"
                exit 1
                ;;
            genePredToGtf)
                echo "genePredToGtf is not installed. You may install using the below command"
                echo "mamba install ucsc-genepredtogtf"
                exit 1
                ;;
        esac
    fi
done


# get options
while getopts "f:g:t:o:p:h" opt; do
    case $opt in
        f) GENOME=$OPTARG;;
        g) GTF=$OPTARG;;
        t) ISFOFORMS=$OPTARG;;
        o) PREFIX=$OPTARG;;
        p) THREADS=$OPTARG;;
        h) usage;;
        *) usage;;
    esac
done

# check if arguments are provided
if [ $# -eq 0 ] || [ -z "$GENOME" ] || [ -z "$GTF" ] || [ -z "$ISFOFORMS" ] || [ -z "$PREFIX" ] || [ -z "$THREADS" ]; then
    usage
fi

# ---- RUN BEGINS HERE ----
start_time=`date +%s`

if [ -f ${PREFIX}.sam.ok ]
then
    print_log "skipping minimap2_mapping"
else
    # minimap2 mapping of isoforms on the genome
    print_log "minimap2_mapping"
    run_cmd "minimap2 -ax splice -C5 --secondary=no -t ${THREADS} ${GENOME} ${ISFOFORMS} | samtools view -F0x900 > ${PREFIX}.sam  2> ${PREFIX}.sam.log"
    touch ${PREFIX}.sam.ok
fi

if [ -f ${PREFIX}.bed.ok ]
then 
    print_log "skipping sam to bed"
else
    # converting sam to bed  
    print_log "sam to bed"
    run_cmd "paftools.js splice2bed ${PREFIX}.sam > ${PREFIX}.bed 2> ${PREFIX}.bed.log"
    touch ${PREFIX}.bed.ok
fi

if [ -f ${PREFIX}.gtf.ok ]
then
    print_log "skipping bed to gtf"
else
    # converting bed to gtf
    print_log "bed to gtf"
    run_cmd "bedToGenePred ${PREFIX}.bed stdout | genePredToGtf -utr file stdin ${PREFIX}.gtf 2> ${PREFIX}.gtf.log"
    touch ${PREFIX}.gtf.ok
fi

if [ -f ${PREFIX}_gffcompare.ok ]
then
    print_log "skipping gffcompare"
else
    # comparing the converted gtf to reference gtf  
    print_log "compare reference to isoform gtf"
    run_cmd "gffcompare -r ${GTF} -o ${PREFIX}_gffcompare ${PREFIX}.gtf 2> ${PREFIX}_gffcompare.log"
    touch ${PREFIX}_gffcompare.ok
fi

if [ -f ${PREFIX}_junceval.stats.ok ]
then
    print_log "skipping summary"
else
    # spliced mapping summary
    print_log "summary"
    run_cmd "paftools.js junceval ${GTF} ${PREFIX}.sam > ${PREFIX}_junceval.stats 2> ${PREFIX}_junceval.stats.log"
    touch ${PREFIX}_junceval.stats.ok
fi

end_time=`date +%s`
run_time=$((end_time-start_time))
print_log "Total run time: $(seconds_to_hhmmss $run_time)"
# ---- RUN ENDS HERE ----

# EOF