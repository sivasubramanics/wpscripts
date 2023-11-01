#!/bin/bash
# script to filter VCF file
# usage: bash filterVCF.sh -i in.vcf.gz -o output_prefix -t threads
# contact: c.s.sivasubramani[at]gmail.com
# date: 18-10-2023

set -e
set -o pipefail
# set -x

MAX_ALLELES=2
MAF=0.05
MIN_QUAL=25
MAX_MISSING=0.1
MAX_HET=0.1
QD=2
THREADS=2

# start clock
start=$(date +%s)

# print log
print_log(){
    echo "[`date +'%Y-%m-%d %H:%M:%S'`] LOG: $1"
}

# run command
run_cmd(){
    echo "[`date +'%Y-%m-%d %H:%M:%S'`] CMD: $1"
    eval $1
}

# convert seconds to hh:mm:ss
seconds_to_hhmmss() {
  local seconds=$1
  local hours=$(printf "%02d" $((seconds / 3600)))
  local minutes=$(printf "%02d" $(((seconds / 60) % 60)))
  local seconds=$(printf "%02d" $((seconds % 60)))
  echo "$hours:$minutes:$seconds"
}

# usage
usage(){
    echo "Usage: $0 [-h] [-i in.vcf.gz] [-o output_prefix] [-t threads]"
    echo "Options:"
    echo "  -i STR    input vcf.gz file"
    echo "  -o STR    output prefix"
    echo "  -t INT    number of threads [2]"
    echo "  -q INT    minimum quality [25]"
    echo "  -m INT    maximum missingness [0.1]"
    echo "  -a INT    maximum number of alleles [2]"
    echo "  -f INT    minimum minor allele frequency [0.05]"
    echo "  -x INT    minimum QD [2]"
    echo "  -y INT    maximum heterozygosity [0.1]"
    echo "  -h        help"
    exit 1
}

# copy file if different
copy(){
    cmp --silent $1 $2 || run_cmd "cp $1 $2"
}

# print inputs
inputs(){
    echo ""
    print_log "========== Input =========="
    print_log "IN_VCF       : ${IN_VCF}"
    print_log "OUT_PREFIX   : ${OUT_PREFIX}"
    print_log "THREADS      : ${THREADS}"
    print_log "MIN_QUAL     : ${MIN_QUAL}"
    print_log "MAX_MISSING  : ${MAX_MISSING}"
    print_log "MAX_ALLELES  : ${MAX_ALLELES}"
    print_log "MAF          : ${MAF}"
    print_log "MAX_HET      : ${MAX_HET}"
    print_log "QD           : ${QD}"
    print_log "==========================="
    echo ""
}

# print outputs
outputs(){
    echo ""
    print_log "========== Output =========="
    print_log "VCF          : ${OUT_PREFIX}.vcf.gz"
    print_log "VCF index    : ${OUT_PREFIX}.vcf.gz.csi"
    print_log "VCF stats    : ${OUT_PREFIX}.vcf.stats"
    print_log "PLINK BIM    : ${OUT_PREFIX}.bim"
    print_log "PLINK FAM    : ${OUT_PREFIX}.fam"
    print_log "PLINK BED    : ${OUT_PREFIX}.bed"
    print_log "PCA          : ${OUT_PREFIX}.cov"
    print_log "============================"
    echo ""
}

# get options
while getopts "i:o:t:q:m:a:f:x:h" opt; do
    case $opt in
        i) IN_VCF=$OPTARG;;
        o) OUT_PREFIX=$OPTARG;;
        t) THREADS=$OPTARG;;
        q) MIN_QUAL=$OPTARG;;
        m) MAX_MISSING=$OPTARG;;
        a) MAX_ALLELES=$OPTARG;;
        f) MAF=$OPTARG;;
        x) QD=$OPTARG;;
        h) usage;;
        *) usage;;
    esac
done

# check if required arguments are provided
if [ $# -eq 0 ] || [ -z "$IN_VCF" ] || [ -z "$OUT_PREFIX" ] || [ -z "$THREADS" ]; then
    usage
fi

# check if required tools are installed
for tool in bcftools plink2 vcf.py; do
    if ! command -v $tool &> /dev/null; then
        case $tool in
            bcftools) 
                echo "Error: $tool is not installed. Install bcftools using the following command:"
                echo "conda install -c bioconda bcftools"
                ;;
            plink2) 
                echo "Error: $tool is not installed. Install plink2 using the following command:"
                echo "conda install -c bioconda plink2"
                ;;
            vcf.py) 
                echo "Error: $tool is not accessible."
                ;;
        esac
        exit 1
    fi
done


TMP_DIR=${OUT_PREFIX}_intermediate_files
TMP_PREFIX=${TMP_DIR}/$(basename ${OUT_PREFIX})

inputs

mkdir -p $TMP_DIR

NUM_SAMPLES=$(bcftools query -l ${IN_VCF} | wc -l)

# This script adds QD and DP tags to the multisample VCF file
QD_VCF=${TMP_PREFIX}.qd.vcf.gz
if [ ! -f ${QD_VCF}.okay ]; then
    run_cmd "vcf.py add_QD -v ${IN_VCF} -o ${QD_VCF} -t ${THREADS} > ${TMP_PREFIX}.qd.vcf.py.log 2>&1"
    touch ${QD_VCF}.okay
else
    print_log "${QD_VCF} already exists. Skipping..."
fi

# Filtering the VCF for mnMAF 0.05, site coverage 0.9, QD 2, minimum quality 25, only biallelic, remove indels 
GTMAX=$(awk '{print ($1-$2)*100}' <<< "1 $MAX_MISSING")
GTHET=$(awk '{print $1*100}' <<< "$MAX_HET")
GTMAF=$(awk '{print $1*100}' <<< "$MAF")
OUTFILE=${TMP_PREFIX}.qd.gt${GTMAX}maf${GTMAF}het${GTHET}al${MAX_ALLELES}qd${QD}.snps.vcf.gz
if [ ! -f ${OUTFILE}.okay ]; then
    run_cmd "bcftools view -i \"MAF >= ${MAF} && F_MISSING <= ${MAX_MISSING} && QD >=${QD} && QUAL >= ${MIN_QUAL} && (COUNT(GT='0/1')/${NUM_SAMPLES}) <= ${MAX_HET}\" --types snps -M ${MAX_ALLELES} --threads ${THREADS} -Oz -o ${OUTFILE} ${TMP_PREFIX}.qd.vcf.gz"
    touch ${OUTFILE}.okay
else
    print_log "${OUTFILE} already exists. Skipping..."
fi

# Indexing
if [ ! -f ${OUTFILE}.csi.okay ]; then
    run_cmd "bcftools index --threads ${THREADS} ${OUTFILE}"
    touch ${OUTFILE}.csi.okay
else
    print_log "${OUTFILE}.csi already exists. Skipping..."
fi

# Basic Stats
STATSFILE=${OUTFILE}.stats
if [ ! -f ${STATSFILE}.okay ]; then
    run_cmd "bcftools stats -s - --threads ${THREADS} ${OUTFILE} > ${STATSFILE}"
    touch ${STATSFILE}.okay
fi

# VCF to Plink
PLINKFILE=${OUTFILE%.vcf.gz}.plk
if [ ! -f ${PLINKFILE}.okay ]; then
    run_cmd "plink2 --vcf ${OUTFILE} --allow-extra-chr --make-bed --out ${PLINKFILE} --vcf-half-call m --threads ${THREADS} > ${PLINKFILE}.log 2>&1"
    touch ${PLINKFILE}.okay
else
    print_log "${PLINKFILE} already exists. Skipping..."
fi

# LD pruning
LDPRUNEFILE=${PLINKFILE}.prune
if [ ! -f ${LDPRUNEFILE}.okay ]; then
    # NUM_SAMPLES=$(awk '{print $2}' ${PLINKFILE}.fam | sort | uniq | wc -l)
    if [ $NUM_SAMPLES -lt 50 ]; then
        print_log "Number of samples is less than 50. Using --bad-ld. Use the data with caution..."
        run_cmd "plink2 --bfile ${PLINKFILE} --indep-pairwise 10'kb' 1 0.5 --out ${LDPRUNEFILE} --allow-extra-chr --threads ${THREADS} --bad-ld > ${LDPRUNEFILE}.log 2>&1"
        run_cmd "plink2 --bfile ${PLINKFILE} --extract ${LDPRUNEFILE}.prune.in --make-bed --out ${LDPRUNEFILE} --allow-extra-chr --threads ${THREADS} > ${LDPRUNEFILE}.log 2>&1"
        touch ${LDPRUNEFILE}.okay
    else
        run_cmd "plink2 --bfile ${PLINKFILE} --indep-pairwise 10'kb' 1 0.5 --out ${LDPRUNEFILE} --allow-extra-chr --threads ${THREADS} > ${LDPRUNEFILE}.log 2>&1"
        run_cmd "plink2 --bfile ${PLINKFILE} --extract ${LDPRUNEFILE}.prune.in --make-bed --out ${LDPRUNEFILE} --allow-extra-chr --threads ${THREADS} > ${LDPRUNEFILE}.log 2>&1"
        touch ${LDPRUNEFILE}.okay
    fi
else
    print_log "${LDPRUNEFILE} already exists. Skipping..."
fi

# Plink to VCF
PLINK2VCF=${LDPRUNEFILE}.vcf.gz
if [ ! -f ${PLINK2VCF}.okay ]; then
    run_cmd "plink2 --bfile ${LDPRUNEFILE} --recode vcf bgz --out ${PLINK2VCF%.vcf.gz} --allow-extra-chr --threads ${THREADS} > ${PLINK2VCF}.log 2>&1"
    touch ${PLINK2VCF}.okay
else
    print_log "${PLINK2VCF} already exists. Skipping..."
fi

# Indexing
if [ ! -f ${PLINK2VCF}.csi.okay ]; then
    run_cmd "bcftools index --threads ${THREADS} ${PLINK2VCF}"
    touch ${PLINK2VCF}.csi.okay
else
    print_log "${PLINK2VCF}.csi already exists. Skipping..."
fi

# Stats
STATSFILE=${PLINK2VCF%.gz}.stats
if [ ! -f ${STATSFILE}.okay ] || [ ! -f ${STATSFILE} ]; then
    run_cmd "bcftools stats -s - --threads ${THREADS} ${PLINK2VCF} > ${STATSFILE}"
    touch ${STATSFILE}.okay
else
    print_log "${STATSFILE} already exists. Skipping..."
fi

# PCA
PCAFILE=${PLINK2VCF%.vcf.gz}.cov
if [ ! -f ${PCAFILE}.okay ]; then
    run_cmd "plink2 --bfile ${LDPRUNEFILE} --pca --out ${PCAFILE} --allow-extra-chr --threads ${THREADS} > ${PCAFILE}.log 2>&1"
    cut -f2,3,4,5 ${PCAFILE}.eigenvec | awk '{if(NR==1){print "taxa\tPC1\tPC2\tPC3"}else{print}}' > ${PCAFILE}
    touch ${PCAFILE}.okay
else
    print_log "${PCAFILE} already exists. Skipping..."
fi

# Copying the results
copy ${PLINK2VCF} ${OUT_PREFIX}.vcf.gz
copy ${PLINK2VCF}.csi ${OUT_PREFIX}.vcf.gz.csi
copy ${STATSFILE} ${OUT_PREFIX}.vcf.stats
copy ${LDPRUNEFILE}.bim ${OUT_PREFIX}.bim
copy ${LDPRUNEFILE}.fam ${OUT_PREFIX}.fam
copy ${LDPRUNEFILE}.bed ${OUT_PREFIX}.bed
copy ${PCAFILE} ${OUT_PREFIX}.cov

outputs

# end clock
end=$(date +%s)
# print runtime
runtime=$((end-start))
print_log "Finished in $(seconds_to_hhmmss $runtime)"

# EOF