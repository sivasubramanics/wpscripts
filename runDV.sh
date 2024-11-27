#!/usr/bin/env bash
# script: runDV.sh
# description: run DeepVariant using Clara Parabricks
# contact: siva.selvanayagam[at]wur.nl
# date: 2023-08-31
# version: 0.0.1

set -u
set -e

sif=/lustre/BIF/nobackup/selva001/work/singularity_imgs/clara-parabricks_4.1.2.sif
threads=2
is_gpu=false
pid=$$

# start clock
start=$(date +%s)

# print log
print_log(){
    echo "[`date +'%Y-%m-%d %H:%M:%S'`] $1"
}

# run command
run_cmd(){
    echo "[`date +'%Y-%m-%d %H:%M:%S'`] CMD: $1"
    eval $1
}

# check if file exists
is_exists(){
    local file=$1
    if [ ! -f $file ]; then
        print_log "ERROR: file not found: $file"
        exit 1
    fi
}

# convert seconds to hh:mm:ss
seconds_to_hhmmss(){
    local seconds=$1
    local hours=$(printf "%02d" $((seconds / 3600)))
    local minutes=$(printf "%02d" $(((seconds / 60) % 60)))
    local seconds=$(printf "%02d" $((seconds % 60)))
    echo "$hours:$minutes:$seconds"
}

# usage
usage(){
    echo "Usage: $0 [-h] [-b bam] [-r ref.fa] [-o out_dir] [-p threads]"
    echo "Options:"
    echo "  -s STR    Clara Parabricks SIF file (eg: $sif)"
    echo "  -b STR    bam/cram file"
    echo "  -r STR    reference genome"
    echo "  -i STR    interval list (bed format)"
    echo "  -o STR    output prefix"
    echo "  -p INT    number of threads"
    echo "  -g        use GPU"
    echo "  -h        help"
    exit 1
}

# get options
while getopts "s:b:r:i:o:p:hg" opt; do
    case $opt in
        s) sif=$OPTARG;;
        b) bam=$OPTARG;;
        r) ref=$OPTARG;;
        i) intrbed=$OPTARG;;
        o) out_prefix=$OPTARG;;
        p) threads=$OPTARG;;
        h) usage;;
        g) is_gpu=true;;
        *) usage;;
    esac
done

if [ $# -eq 0 ]; then
    usage
fi

# check options
if [[ -z $sif || -z $bam || -z $ref || -z $out_prefix || -z $threads ]]; then
    usage
fi

# check if files exist
is_exists $sif
is_exists $bam
is_exists $ref

# check tools
for tool in singularity; do
    if ! command -v $tool &> /dev/null; then
        print_log "ERROR: $tool is not installed"
        exit 1
    fi
done


out_dir=$(dirname $out_prefix)
out_tmp="$out_dir/tmp_$pid"
# create output directory
mkdir -p $out_dir
mkdir -p $out_tmp

sif=$(realpath $sif)
in_bam_path=$(dirname $(realpath $bam))
in_bam_name=$(basename $bam)
in_ref_path=$(dirname $(realpath $ref))
in_ref_name=$(basename $ref)
# if interval bed file is provided
if [ ! -z $intrbed ]; then
    is_exists $intrbed
    in_intrbed_path=$(dirname $(realpath $intrbed))
    in_intrbed_name=$(basename $intrbed)
fi
out_vcf_dir=$(realpath $out_dir/)
tmp_dir=$(realpath $out_tmp/)
out_gvcf_name=$(basename $out_prefix).clara.g.vcf.gz
out_vcf_name=$(basename $out_prefix).clara.vcf
log_file=$(basename $out_prefix).clara.log
run_log=$(basename $out_prefix).clara.run.log

dv_cmd="singularity run"
if $is_gpu; then
    dv_cmd="$dv_cmd --nv"
fi
dv_cmd="$dv_cmd -B $in_bam_path:/bams"
dv_cmd="$dv_cmd -B $in_ref_path:/genome"
if [ ! -z $intrbed ]; then
    dv_cmd="$dv_cmd -B $in_intrbed_path:/intervals"
    dv_cmd="$dv_cmd --interval-file /intervals/$in_intrbed_name"
fi
dv_cmd="$dv_cmd -B $out_vcf_dir:/vcfs"
dv_cmd="$dv_cmd -B $tmp_dir:/tmp"
dv_cmd="$dv_cmd $sif"
dv_cmd="$dv_cmd pbrun deepvariant"
dv_cmd="$dv_cmd --ref /genome/$in_ref_name"
dv_cmd="$dv_cmd --in-bam /bams/$in_bam_name"
dv_cmd="$dv_cmd --out-variants /vcfs/$out_gvcf_name"
# dv_cmd="$dv_cmd --pb-model-file /usr/local/parabricks/binaries/model/75/shortread/deepvariant.eng"
dv_cmd="$dv_cmd --gvcf"
dv_cmd="$dv_cmd --logfile /vcfs/$log_file"
# dv_cmd="$dv_cmd --num-cpu-threads-per-stream $threads"
dv_cmd="$dv_cmd --num-gpus 2"

run_cmd "$dv_cmd"

print_log "compressing vcf"
run_cmd "bgzip --threads $threads -f $out_vcf_dir/$out_vcf_name"

# remove tmp files
run_cmd "rm -rf $out_tmp"

# stop clock
end=$(date +%s)

# print log
print_log "elapsed time: $(seconds_to_hhmmss $((end - start)))"

# END