#!/bin/bash

set -u
set -e

start=$(date +%s)

print_log(){
    echo "[`date +'%Y-%m-%d %H:%M:%S'`] $1"
}

run_cmd(){
    echo "[`date +'%Y-%m-%d %H:%M:%S'`] CMD: $1"
    eval $1
}

usage(){
    echo "Usage: $0 [-h] [-f fq.list] [-r ref.fa] [-o output] [-p threads] [-k kmer_size] [-w window_size]"
    echo "Options:"
    echo "  -f STR    list of fastq files"
    echo "  -r STR    reference genome"
    echo "  -o STR    output file"
    echo "  -p INT    number of threads"
    echo "  -k INT    kmer size"
    echo "  -w INT    window size"
    echo "  -h        help"
    exit 1
}

check_exec(){
    if ! command -v $1 &> /dev/null; then
        echo "Error: $1 is not installed"
        exit 1
    fi
}

while getopts "f:r:o:p:k:w:h" opt; do
    case $opt in
        f) fq_list=$OPTARG;;
        r) ref=$OPTARG;;
        o) output=$OPTARG;;
        p) threads=$OPTARG;;
        k) kmer_size=$OPTARG;;
        w) window_size=$OPTARG;;
        h) usage;;
        *) usage;;
    esac
done

if [ -z "$fq_list" ] || [ -z "$ref" ] || [ -z "$output" ] || [ -z "$threads" ] || [ -z "$kmer_size" ] || [ -z "$window_size" ]; then
    usage
fi

print_log "Checking tools"
check_exec kmc 
check_exec $TOOLS_PATH/kmersGWAS/bin/kmers_add_strand_information

print_log "Creating output directory $output"
mkdir -p $output

print_log "KMC canonized counting"
run_cmd "kmc -k$kmer_size -t$threads -ci0 @$fq_list $output/canon $output" > $output/kmc.canon.log 2>&1

print_log "KMC no canonized counting (all k-mers)"
run_cmd "kmc -k$kmer_size -t$threads -ci0 -b @$fq_list $output/all $output" > $output/kmc.all.log 2>&1

print_log "Adding strand information"
run_cmd "$TOOLS_PATH/kmersGWAS/bin/kmers_add_strand_information -c $output/canon -n $output/all -k $kmer_size -o $output/kmc31" > $output/kmers_add_strand_information.log 2>&1

print_log "IBScpp variant calling"
run_cmd "$TOOLS_PATH/IBSpy/IBScpp/build/IBScpp -d $output/kmc31 -r $ref -p $threads -k $kmer_size -w $window_size > $output/variations.tsv" > $output/IBScpp.variations.log 2>&1

end=$(date +%s)
runtime=$((end-start))

print_log "Runtime: $runtime seconds"

