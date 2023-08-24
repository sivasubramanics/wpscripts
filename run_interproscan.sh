#!/bin/bash
# script: run_interproscan.sh
# description: run interproscan
# contact: siva.selvanayagam[at]wur.nl
# date: 2023-08-21
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
    echo "Usage: $0 [-h] [-f in.fasta] [-o output_dir] [-t threads]"
    echo "Options:"
    echo "  -f STR    input fasta file"
    echo "  -o STR    output directory"
    echo "  -t INT    number of threads"
    echo "  -h        help"
    exit 1
}

# get options
while getopts "f:o:t:h" opt; do
    case $opt in
        f) in_fasta=$OPTARG;;
        o) out_dir=$OPTARG;;
        t) THREADS=$OPTARG;;
        h) usage;;
        *) usage;;
    esac
done

# check if required tools are installed
for tool in $TOOLS_PATH/interproscan-5.59-91.0/interproscan.sh faSplit; do
    if ! command -v $tool &> /dev/null; then
        echo "ERROR: $tool not found"
        exit 1
    fi
done

# check if input file exists
if [ ! -f "$in_fasta" ]; then
    echo "ERROR: $in_fasta file not found"
    exit 1
fi

# check if output directory exists
if [ -d "$out_dir" ]; then
    echo "ERROR: $out_dir directory already exists"
    exit 1
fi

# rm -rf $out_dir
# create output directory
mkdir -p $out_dir

# split fasta file
print_log "Splitting fasta file"
# split fasta file to 1000 sequences per file
awk '$1~/^>/{print substr($1,2);}' $in_fasta > $in_fasta.ids
split -l 100 $in_fasta.ids -d $out_dir/split
ls $out_dir/split* > $(basename $out_dir).list
for i in $(cat $(basename $out_dir).list); do
    echo "faSomeRecords $in_fasta $i $i.fasta"
done | parallel -j $THREADS

# run interproscan
print_log "Running interproscan"
for i in $(cat $(basename $out_dir).list); do
    tmp_dir=$i.tmp
    mkdir -p $tmp_dir
    ipr_cmd="$TOOLS_PATH/interproscan-5.59-91.0/interproscan.sh \
            -appl TIGRFAM,SUPERFAMILY,PANTHER,Gene3D,Coils,Pfam,MobiDBLite \
            -iprlookup \
            -goterms \
            -pa \
            -i $i.fasta \
            -d $tmp_dir \
            -cpu 4 &> $i.iprscan.log 2>&1"
    echo $ipr_cmd >> $out_dir/ipr_cmd.sh
done
parallel -j $THREADS < $out_dir/ipr_cmd.sh

# merge results
print_log "Merging results"
rm -f $(basename $out_dir)_iprscan.ann.tsv
rm -f $(basename $out_dir)_iprscan.gff3
for i in $(cat $(basename $out_dir).list); do
    tmp_dir=$i.tmp
    if [ ! -f $(basename $out_dir)_iprscan.ann.tsv ]; then
        cat $tmp_dir/*.tsv >> $(basename $out_dir)_iprscan.ann.tsv
        cat $tmp_dir/*.gff3 >> $(basename $out_dir)_iprscan.gff3
    else
        cat $tmp_dir/*.tsv >> $(basename $out_dir)_iprscan.ann.tsv
        tail -n +3 $tmp_dir/*.gff3 >> $(basename $out_dir)_iprscan.gff3
    fi
done

# remove temporary files
print_log "Removing temporary files"
rm -rf $out_dir

# stop clock
end=$(date +%s)

# print log
print_log "Total time elapsed: $(seconds_to_hhmmss $((end-start)))"

# end