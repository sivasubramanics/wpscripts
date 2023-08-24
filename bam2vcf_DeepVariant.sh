#!/usr/bin/env bash
# script: bam2vcf_DeepVariant.sh
# description: call variations from bam file using DeepVariant
# contact: siva.selvanayagam[at]wur.nl
# date: 2023-07-18
# version: 0.0.1

set -u
set -e

sif="/lustre/BIF/nobackup/selva001/work/singularity_imgs/dv_1.5.0.sif"
BIN_VERSION="1.5.0"
threads=1
chrom="all"

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
    echo "  -s STR    deepvariant SIF file (eg: $sif)"
    echo "  -b STR    bam file"
    echo "  -r STR    reference genome"
    echo "  -o STR    output directory"
    echo "  -c STR    chromosome name [or] file containing chromosome names [or] comma separated chromosome names"
    echo "  -p INT    number of threads"
    echo "  -h        help"
    exit 1
}

# get options
while getopts "s:b:r:o:p:c:h" opt; do
    case $opt in
        s) sif=$OPTARG;;
        b) bam=$OPTARG;;
        r) ref=$OPTARG;;
        o) out_dir=$OPTARG;;
        p) threads=$OPTARG;;
        c) chrom=$OPTARG;;
        h) usage;;
        *) usage;;
    esac
done

# check options
if [ $# -eq 0 ] || [ -z "$bam" ] || [ -z "$ref" ] || [ -z "$out_dir" ]; then
    usage
fi

# check tools
for tool in samtools bcftools singularity parallel; do
    if ! command -v $tool &> /dev/null; then
        echo "Error: $tool is not installed"
        exit 1
    fi
done

# check input files
if [ ! -f $sif ]; then
    echo "Error: $sif not found"
    exit 1
fi

# check input files
if [ ! -f $ref.fai ]; then
    print_log "indexing reference genome"
    run_cmd "samtools faidx $ref"
fi

if [ ! -f $bam.bai ]; then
    print_log "indexing bam file"
    run_cmd "samtools index $bam"
fi

# create output directory
mkdir -p $out_dir
mkdir -p $out_dir/genome
mkdir -p $out_dir/logs
mkdir -p $out_dir/vcfs
mkdir -p $out_dir/out
mkdir -p $out_dir/sh

if [ "$chrom" == "all" ]; then
    chroms=$(cut -f1 $ref.fai | tr '\n' ' ')
elif [ -f $chrom ]; then
    chroms=$(cut -f1 $chrom | tr '\n' ' ')
else
    chroms=$(echo $chrom | tr ',' ' ')
fi

print_log "remove old files"
run_cmd "rm -f $out_dir/sh/dv.sh $out_dir/vcfs.list $out_dir/gvcfs.list"

sif=$(realpath $sif)
in_bam_path=$(dirname $(realpath $bam))
in_bam_name=$(basename $bam)
in_ref_path=$(dirname $(realpath $ref))
in_ref_name=$(basename $ref)
out_vcf_dir=$(realpath $out_dir/vcfs)

print_log "generate non-gap regions for reference genome"
run_cmd "seqtk cutN $ref | grep \">\" | sed 's/>//' > $out_dir/genome/all.list"

for chr in $chroms; do
    reg_count=0
    
    print_log "creating output directory for chromosome: $chr"
    run_cmd "mkdir -p $out_dir/vcfs/$chr"
    run_cmd "mkdir -p $out_dir/logs/dv_$chr"
    run_cmd "rm -f $out_dir/vcfs/$chr.vcfs.list $out_dir/vcfs/$chr.gvcfs.list"

    print_log "creating list of regions for chromosome: $chr"
    run_cmd "grep \"$chr:\" $out_dir/genome/all.list > $out_dir/genome/$chr.list"

    for reg in $(cat $out_dir/genome/$chr.list); do
        reg_count=$((reg_count+1))
        run_cmd "mkdir -p $out_dir/tmp/$chr/$reg_count"
        in_tmp_dir=$(realpath $out_dir/tmp/$chr/$reg_count)
        out_vcf_name=$chr/$reg_count.vcf.gz
        out_gvcf_name=$chr/$reg_count.g.vcf.gz

        dv_cmd="singularity run -B /usr/lib/locale/:/usr/lib/locale/"
        dv_cmd="${dv_cmd} -B $in_tmp_dir:/tmp" 
        dv_cmd="${dv_cmd} -B $in_bam_path:/bam"
        dv_cmd="${dv_cmd} -B $in_ref_path:/ref"
        dv_cmd="${dv_cmd} -B $out_vcf_dir:/vcf"
        dv_cmd="${dv_cmd} $sif /opt/deepvariant/bin/run_deepvariant"
        dv_cmd="${dv_cmd} --model_type=WGS"
        dv_cmd="${dv_cmd} --ref=/ref/$in_ref_name"
        dv_cmd="${dv_cmd} --reads=/bam/$in_bam_name"
        dv_cmd="${dv_cmd} --output_vcf=/vcf/$out_vcf_name"
        dv_cmd="${dv_cmd} --output_gvcf=/vcf/$out_gvcf_name"
        dv_cmd="${dv_cmd} --intermediate_results_dir /tmp"
        dv_cmd="${dv_cmd} --regions $reg"
        dv_cmd="${dv_cmd} --num_shards=4"
        dv_cmd="${dv_cmd} --call_variants_extra_args batch_size=5120"
        dv_cmd="${dv_cmd} > $out_dir/logs/dv_$chr/$reg_count.log 2>&1"

        echo $dv_cmd >> $out_dir/sh/dv.sh
        echo $out_vcf_dir/$out_vcf_name >> $out_dir/vcfs/$chr.vcfs.list
        echo $out_vcf_dir/$out_gvcf_name >> $out_dir/vcfs/$chr.gvcfs.list
    done
done

print_log "running deepvariant in parallel"
parallel -j $((threads / 2 )) < $out_dir/sh/dv.sh > $out_dir/logs/dv.log 2>&1

for chr in $chroms; do
    print_log "merging vcf files: $chr"
    run_cmd "bcftools concat --threads $threads -f $out_dir/vcfs/$chr.vcfs.list -Oz -o $out_dir/out/$chr.vcf.gz > $out_dir/logs/vcf_concat.log 2>&1"
    run_cmd "bcftools index  --threads $threads -f -t $out_dir/out/$chr.vcf.gz"
    run_cmd "bcftools concat --threads $threads -f $out_dir/vcfs/$chr.gvcfs.list -Oz -o $out_dir/out/$chr.g.vcf.gz > $out_dir/logs/gvcf_concat.log 2>&1"
    run_cmd "bcftools index --threads $threads -f -t $out_dir/out/$chr.g.vcf.gz"
done

print_log "removing tmp files"
run_cmd "rm -rf $out_dir/tmp $out_dir/genome "


# end clock
end=$(date +%s)
runtime=$((end-start))

# print runtime
print_log "Finished in $(seconds_to_hhmmss $runtime)"

