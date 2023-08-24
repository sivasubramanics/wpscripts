#!/bin/bash
# script: run_fastp.sh
# description: run fastp on all fastq files in a directory
# contact: siva.selvanayagam[at]wur.nl
# date: 2023-07-21
# version: 0.1

set -u
set -e

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
    echo "Usage: $0 [-h] [-i in1.fastq] -j [in2.fastq] [-o output_dir] [-t threads]"
    echo "Options:"
    echo "  -i STR    input fastq file read1"
    echo "  -j STR    input fastq file read2"
    echo "  -o STR    output directory"
    echo "  -t INT    number of threads"
    echo "  -h        help"
    exit 1
}

# get options
while getopts "i:j:o:t:h" opt; do
    case $opt in
        i) in_fastq_one=$OPTARG;;
        j) in_fastq_two=$OPTARG;;
        o) out_dir=$OPTARG;;
        t) THREADS=$OPTARG;;
        h) usage;;
        *) usage;;
    esac
done

# check if required tools are installed
for tool in fastp; do
    if ! command -v $tool &> /dev/null; then
        echo "Error: $tool is not installed"
        exit 1
    fi
done

# check if required arguments are provided
if [ -z "$in_fastq_one" ]; then
    echo "Error: input fastq read 1 file not provided"
    usage
fi
if [ -z "$in_fastq_two" ]; then
    echo "Error: input fastq read 2 file not provided"
    usage
fi
if [ -z "$out_dir" ]; then
    echo "Error: output directory not provided"
    usage
fi
if [ -z "$THREADS" ]; then
    print_log "number of threads not provided. Using 1 thread"
    THREADS=1
fi

# check if output directory exists
if [ ! -d $out_dir ]; then
    mkdir -p $out_dir
fi

if [ ! -f $in_fastq_one ]; then
    echo "Error: $in_fastq_one does not exist"
    exit 1
fi
if [ ! -f $in_fastq_two ]; then
    echo "Error: $in_fastq_two does not exist"
    exit 1
fi

base_one=$(basename $in_fastq_one)
base_two=$(basename $in_fastq_two)
fwd_fq=$out_dir/${base_one/.fastq.gz/.flt.fastq.gz}
fwd_fq=${fwd_fq/.fq.gz/.flt.fq.gz}
rev_fq=$out_dir/${base_two/.fastq.gz/.flt.fastq.gz}
rev_fq=${rev_fq/.fq.gz/.flt.fq.gz}

html=$out_dir/fastp_out.html
json=$out_dir/fastp_out.json
run_cmd "fastp --thread $THREADS \
        --detect_adapter_for_pe \
        --correction \
        --qualified_quality_phred 20 \
        --unqualified_percent_limit 30 \
        --cut_tail \
        --cut_tail_window_size 1 \
        --cut_tail_mean_quality 3 \
        --cut_front \
        --cut_front_window_size 1 \
        --cut_front_mean_quality 3 \
        --length_required 35 \
        --html $html \
        --json $json \
        --in1 $in_fastq_one --in2 $in_fastq_two \
        --out1 $fwd_fq --out2 $rev_fq &> $out_dir/fastp.log 2>&1"


# end clock
end=$(date +%s)
runtime=$((end-start))

# print runtime
echo "Runtime: $(seconds_to_hhmmss $runtime)"


