#!/usr/bin/env bash

set -u
set -e

THREADS=2
RESUME=false
IPRSCAN_EXE=/lustre/BIF/nobackup/selva001/work/tools/interproscan-5.64-96.0/interproscan.sh

print_log(){
    echo "[`date +'%Y-%m-%d %H:%M:%S'`] $1"
}

run_cmd(){
    echo "[`date +'%Y-%m-%d %H:%M:%S'`] CMD: $1"
    eval $1
}

is_exists(){
    local file=$1
    if [ ! -f $file ]; then
        print_log "ERROR: file not found: $file"
        exit 1
    fi
}

seconds_to_hhmmss(){
    local seconds=$1
    local hours=$(printf "%02d" $((seconds / 3600)))
    local minutes=$(printf "%02d" $(((seconds / 60) % 60)))
    local seconds=$(printf "%02d" $((seconds % 60)))
    echo "$hours:$minutes:$seconds"
}

usage(){
    echo "Usage: $0 [-h] [-t isoforms.fna] [-o out_dir] [-p threads] [-r]"
    echo "Options:"
    echo "  -t STR    isoforms/transcripts fasta file"
    echo "  -o STR    output directory"
    echo "  -p INT    number of threads"
    echo "  -r        resume pipeline"
    echo "  -h        help"
    exit 1
}

while getopts "t:o:p:a:m:hr" opt; do
    case $opt in
        t) TRANSCRIPTOME=$OPTARG;;
        o) OUT_DIR=$OPTARG;;
        p) THREADS=$OPTARG;;
        a) PROTEOME=$OPTARG;;
        m) TR2AAMAP=$OPTARG;;
        r) RESUME=true;;
        h) usage;;
        *) usage;;
    esac
done

if [ $# -eq 0 ]; then
    usage
fi

# check if transcriptome file exists
is_exists $TRANSCRIPTOME

# check if output directory exists
if [ ! -d $OUT_DIR ]; then
    run_cmd "mkdir -p $OUT_DIR"
fi

# check tools
for tool in $IPRSCAN_EXE TransDecoder.LongOrfs TransDecoder.Predict; do
    if ! command -v $tool &> /dev/null; then
        print_log "ERROR: $tool is not installed"
        exit 1
    fi
done

remove_duplicates(){
    print_log "Removing duplicate sequences..."
    run_cmd "seqkit rmdup -s -i $TRANSCRIPTOME -o $OUT_DIR/transcripts.fna"
    run_cmd "touch $OUT_DIR/transcripts.fna.ok"
}

if $RESUME; then
    if [ ! -f $OUT_DIR/transcripts.fna.ok ]; then
        remove_duplicates
    else
        print_log "skipping duplicate removal."
    fi
else
    remove_duplicates
fi


run_annotation(){
    print_log "Running interpro annotation..."
    run_cmd "$IPRSCAN_EXE -i $TRANSCRIPTOME -t rna -f tsv -o $OUT_DIR/interproscan.tsv -dp -pa -goterms -iprlookup -cpu $THREADS -appl Pfam,PANTHER"
    run_cmd "touch $OUT_DIR/interproscan.tsv.ok"
}

if $RESUME; then
    if [ ! -f $OUT_DIR/interproscan.tsv.ok ]; then
        run_annotation
    else
        print_log "skipping interpro annotation."
    fi
else
    run_annotation
fi

run_eggnog(){
    print_log "Running eggnog annotation..."
    mkdir -p $OUT_DIR/eggnog/tmp
    run_cmd "emapper.py -i $TRANSCRIPTOME --output $OUT_DIR/eggnog --cpu $THREADS --data_dir /lustre/BIF/nobackup/selva001/work/tools/eggnog-mapper-2.1.2/data --override --seed_ortholog_evalue 1e-5 --seed_ortholog_score 40 --target_ortholog_evalue 1e-5 --target_ortholog_score 40 --go_evalue 1e-5 --go_weight 5 --no_file_comments"
}

tr2aa(){
    print_log "Translating transcripts to proteins..."
    run_cmd "seqkit translate -i $OUT_DIR/transcripts.fna -o $OUT_DIR/proteins.faa"
    run_cmd "touch $OUT_DIR/proteins.faa.ok"
}