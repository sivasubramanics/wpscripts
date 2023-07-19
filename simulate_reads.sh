#!/bin/bash
# script: simulate_reads.sh
# description: simulate reads from reference genome
# contact: siva.selvanayagam[at]wur.nl
# date: 2023-07-19
# version: 0.0.1

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

get_cov(){
    local mean=$1
    local std_dev=${2:-0.3}
    local min=$(echo "$mean - ($mean * $std_dev)" | bc)
    local max=$(echo "$mean + ($mean * $std_dev)" | bc)
    local range=$(echo "$max - $min + 1" | bc)
    local rand=$(awk 'BEGIN{srand(); print rand()}')
    echo $(echo "$min + $rand * ($max - $min)" | bc)
}

# function to convert seconds to human readable time
seconds_to_hhmmss() {
  local seconds=$1
  local hours=$(printf "%02d" $((seconds / 3600)))
  local minutes=$(printf "%02d" $(((seconds / 60) % 60)))
  local seconds=$(printf "%02d" $((seconds % 60)))
  echo "$hours:$minutes:$seconds"
}
# usage message
usage(){
    echo "Usage: $0 [-h] [-t type] [-r ref] [-o out_dir] [-p threads]"
    echo "Options:"
    echo "  -t STR    type of reads (WGS or RNA)"
    echo "  -r STR    reference genome name (one of: 'yeast', 'fruitfly', 'human', 'mouse', 'arabidopsis', 'rice', 'zebrafish', 'maize', 'tomato', 'lettuce')"
    echo "  -o STR    output directory"
    echo "  -n INT    number of samples"
    echo "  -c FLT    coverage"
    echo "  -h        help"
    exit 1
}

n_samples=1
n_reps=3

# parse arguments
while getopts "t:r:o:n:c:h" opt; do
    case $opt in
        t) type=$OPTARG;;
        r) ref_name=$OPTARG;;
        o) out_dir=$OPTARG;;
        n) n_samples=$OPTARG;;
        c) coverage=$OPTARG;;
        h) usage;;
        *) usage;;
    esac
done

# check arguments
if [ $# -eq 0 ] || [ -z "$type" ] || [ -z "$ref_name" ] || [ -z "$out_dir" ] || [ -z "$coverage" ]; then
    usage
fi

# check tools
for tool in datasets randomreads.sh; do
    if ! command -v $tool &> /dev/null; then
        case $tool in
            datasets)
                echo "Error: datasets is not installed. Please install datasets using the following command:"
                echo "conda install -c conda-forge ncbi-datasets-cli"
                ;;
            randomreads.sh)
                echo "Error: randomreads.sh is not installed. Please install randomreads.sh using the following command:"
                echo "conda install -c bioconda bbmap"
                ;;
        esac
        exit 1
    fi
done

case $ref_name in
    yeast)
        ref="GCF_000146045.2"
        ;;
    fruitfly)
        ref="GCF_000001215.4"
        ;;
    human)
        ref="GCF_000001405.40"
        ;;
    mouse)
        ref="GCF_000001635.27"
        ;;
    arabidopsis)
        ref="GCF_000001735.4"
        ;;
    rice)
        ref="GCF_001433935.1"
        ;;
    zebrafish)
        ref="GCF_000002035.6"
        ;;
    maize)
        ref="GCF_902167145.1"
        ;;
    tomato)
        ref="GCF_000188115.4"
        ;;
    lettuce)
        ref="GCF_002870075.4"
        ;;
    *)
        echo "Error: reference genome $ref_name is not supported. Please choose one from 'yeast', 'fruitfly', 'human', 'mouse', 'arabidopsis', 'rice', 'zebrafish', 'maize', 'tomato', 'lettuce'"; 
        exit 1
        ;;
esac

# create output directory
{
    mkdir -p $out_dir
    mkdir -p $out_dir/download
    mkdir -p $out_dir/logs
    mkdir -p $out_dir/genome
    mkdir -p $out_dir/reads
    ref_genome_fa="$out_dir/genome/$ref_name.genome.fa"
    ref_rna_fa="$out_dir/genome/$ref_name.rna.fa"
    ref_gff="$out_dir/genome/$ref_name.genome.gff"
}

# download reference genome
if [ ! -f $ref_genome_fa ] || [ ! -f $ref_rna_fa ] || [ ! -f $ref_gff ]; then
    print_log "Downloading reference genome $ref_name"
    run_cmd "datasets download genome accession GCF_000146045.2 --include gff3,rna,genome --filename $out_dir/download/ncbi_dataset.zip --no-progressbar > $out_dir/logs/download.log 2>&1"
    print_log "Unzipping reference genome"
    run_cmd "unzip -o -q $out_dir/download/ncbi_dataset.zip -d $out_dir/download > $out_dir/logs/unzip.log 2>&1"
    print_log "Moving reference genome to $out_dir/genome"
    run_cmd "mv $out_dir/download/ncbi_dataset/data/GCF_000146045.2/$ref*_genomic.fna $ref_genome_fa"
    run_cmd "mv $out_dir/download/ncbi_dataset/data/GCF_000146045.2/genomic.gff $ref_gff"
    run_cmd "mv $out_dir/download/ncbi_dataset/data/GCF_000146045.2/rna.fna $ref_rna_fa"
fi

if [ "$type" == "WGS" ]; then
    ref_fa=$ref_genome_fa
fi
if [ "$type" == "RNA" ]; then
    ref_fa=$ref_rna_fa
fi

for i in $(seq 1 $n_samples); do
    if [ "$type" == "WGS" ]; then
        print_log "Simulating reads for sample $i"
        prefix=$ref_name"_"$type"_"sample_$i
        cov=$(get_cov $coverage 0.1)
        run_cmd "randomreads.sh ref=$ref_fa out1=$out_dir/reads/$prefix.1.fq.gz out2=$out_dir/reads/$prefix.2.fq.gz overwrite=true illuminanames=true paired=true mininsert=100 maxinsert=1000 minlength=100 maxlength=100 coverage=$cov > $out_dir/logs/randomreads_$prefix.log 2>&1"
    fi
    if [ "$type" == "RNA" ]; then
        cov=$(get_cov $coverage 0.5)
        for j in $(seq 1 $n_reps); do
            print_log "Simulating reads for sample $i rep $j"
            prefix=$ref_name"_"$type"_"sample_$i"_"rep_$j
            cov=$(get_cov $cov 0.001)
            run_cmd "randomreads.sh ref=$ref_fa out1=$out_dir/reads/$prefix.1.fq.gz out2=$out_dir/reads/$prefix.2.fq.gz overwrite=true illuminanames=true paired=true mininsert=100 maxinsert=1000 minlength=100 maxlength=100 coverage=$cov > $out_dir/logs/randomreads_$prefix.log 2>&1"
        done
    fi
done

run_cmd "rm -rf ref $out_dir/download"
end=$(date +%s)
runtime=$((end-start))

print_log "Finished in $(seconds_to_hhmmss $runtime)"

