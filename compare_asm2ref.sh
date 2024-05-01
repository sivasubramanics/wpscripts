#!/bin/bash
# This script will do the following steps,
#   1. map the input transcirpt sequences to the reference genome using minimap2
#   2. convert the output to PSL
#   3. filter the PSL file to remove the low quality alignments, and convert to GTF
#   4. compare the GTF to the reference annotation using gffcompare
# USAGE: ./run_alignment.sh -g <genome.fa> -t <transcripts.fa> -a <annotation.gtf> -o <output_dir>
# REQUIREMENTS:
#   1. minimap2
#   2. gffcompare
#   3. samtools
#   4. psl-filter
#   5. sam2psl.py

set -e
set -u

# default number of threads
threads=2
force=0
cutoff_score=200


# function to print log
function prtlog(){
    echo "[`date +'%Y-%m-%d %H:%M:%S'`] LOG: $1"
}

# function to run command
function run_cmd(){
    cmd=$1
    prtlog "Running command: $cmd"
    eval $cmd
}

# function to print usage
function usage(){
    echo "USAGE: $0 -g <genome.fa> -t <transcripts.fa> -a <annotation.gtf> -o <output_dir> [-p <num_threads>] [-f]"
    echo "OPTIONS:"
    echo "  -g  <genome.fa>        : reference genome in fasta format"
    echo "  -t  <transcripts.fa>   : input transcripts in fasta format"
    echo "  -a  <annotation.gtf>   : reference annotation in GTF format"
    echo "  -o  <output_dir>       : output directory"
    echo "  -p  <num_threads>      : number of threads to use (default: 2)"
    echo "  -f                     : force re-run the pipeline even if the output files exist"
    echo "  -s  <cutoff_score>     : minimum PSL/GTF score to keep 0-1000 (default: 200)"
    echo "  -h                     : print this help message"
    exit 1
}

# function to convert seconds to hh:mm:ss
function sec2hhmmss(){
    local T=$1
    local H=$((T/3600))
    local M=$((T%3600/60))
    local S=$((T%60))
    printf "%02d:%02d:%02d\n" $H $M $S
}

for tool in minimap2 gffcompare samtools psl-filter sam2psl.py; do
    if ! command -v $tool &> /dev/null; then
        echo "ERROR: $tool is not in the PATH"
        exit 1
    fi
done

# parse arguments
while getopts "g:t:a:o:hfp:s:" opt; do
    case $opt in
        g) genome=$OPTARG;;
        t) transcripts=$OPTARG;;
        a) annotation=$OPTARG;;
        o) outdir=$OPTARG;;
        p) threads=$OPTARG;;
        f) force=1;;
        h) usage;;
        \?) usage;;
    esac
done

# check if required arguments are provided
if [ $# -eq 0 ]; then
    usage
fi

if [ -z $genome ] || [ -z $transcripts ] || [ -z $annotation ] || [ -z $outdir ]; then
    echo "ERROR: missing required arguments"
    usage
fi

# ----- RUN BIGINS HERE -----
start_time=`date +%s`

# if force is set then delete the output directory
if [ $force -eq 1 ]; then
    prtlog "Removing output directory $outdir"
    run_cmd "rm -rf $outdir"
fi

# create output directory if not exists
run_cmd "mkdir -p $outdir"

sam_file=$outdir/aln.sam
sam_file_chk=$outdir/aln.sam.ok
psl_file=$outdir/aln.psl
psl_file_chk=$outdir/aln.psl.ok
gtf_file=$outdir/aln.flt.gtf
gtf_file_chk=$outdir/aln.flt.gtf.ok
cmp_file=$outdir/cmp.stats
cmp_file_chk=$outdir/cmp.stats.ok
ids_all=$outdir/ids.all.list
ids_mapped=$outdir/ids.mapped.list
ids_unmapped=$outdir/ids.unmapped.list
fa_mapped=$outdir/mapped.fasta
fa_mapped_chk=$outdir/mapped.fasta.ok
fa_unmapped=$outdir/unmapped.fasta
fa_unmapped_chk=$outdir/unmapped.fasta.ok

# map the input transcripts to the reference genome
if [ -f $sam_file_chk ]
then
    prtlog "Skipping mapping, $sam_file already exists"
else
    run_cmd "minimap2 -ax splice -C5 --secondary=no -t $threads $genome $transcripts | samtools view -F0x900 > $sam_file  2> $sam_file.log"
    run_cmd "touch $sam_file_chk"
fi

# convert the sam file to psl
if [ -f $psl_file_chk ]
then
    prtlog "Skipping conversion, $psl_file already exists"
else
    run_cmd "sam2psl.py -i $sam_file -o $psl_file"
    run_cmd "touch $psl_file_chk"
fi

# filter the psl file
if [ -f $gtf_file_chk ]
then
    prtlog "Skipping filtering, $gtf_file already exists"
else
    tmp_out_prefix=${gtf_file%.gtf}
    run_cmd "psl-filter genome_psl -i $psl_file -o $tmp_out_prefix -s $cutoff_score"
    run_cmd "touch $gtf_file_chk"
fi

# compare the filtered GTF to the reference annotation
if [ -f $cmp_file_chk ]
then
    prtlog "Skipping comparison, $cmp_file already exists"
else
    tmp_out_prefix=${cmp_file%.stats}
    run_cmd "gffcompare -r $annotation -G -o $tmp_out_prefix $gtf_file 2> $tmp_out_prefix.log"
    run_cmd "touch $cmp_file_chk"
fi

if [ -f $fa_mapped_chk ]
then
    prtlog "Skipping extraction, $fa_mapped already exists"
else
    run_cmd "awk '{if($0~/^>/){s=substr($1,2); print s}}' $transcripts > $ids_all"
    run_cmd "cut -f1 ${gtf_file%.gtf}.loci | sed '1d' | sort | uniq > $ids_mapped"
    run_cmd "faSomeRecords $transcripts $ids_mapped $fa_mapped"
    run_cmd "touch $fa_mapped_chk"
fi

if [ -f $fa_unmapped_chk ]
then
    prtlog "Skipping extraction, $fa_unmapped already exists"
else
    run_cmd "cat $ids_all $ids_mapped | sort | uniq -u > $ids_unmapped"
    run_cmd "faSomeRecords $transcripts $ids_unmapped -exclude $fa_unmapped"
    run_cmd "touch $fa_unmapped_chk"
fi

# ----- RUN ENDS HERE -----

# print the runtime
end_time=`date +%s`
runtime=$((end_time-start_time))
prtlog "Total runtime: $(sec2hhmmss $runtime)"

# EOF
