#!/usr/bin/env bash
# Run GLnexus on gVCF files

set -u
set -e

threads=2
pid=$$
GL_BIN="/lustre/BIF/nobackup/selva001/work/tools/bin/glnexus_cli"
tmp_dir=""

start_time=$(date +%s)

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

num_lines(){
    local file=$1
    wc -l $file | cut -f1 -d' '
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
    echo "Usage: $0 [-h] [-o output.vcf.gz] [-t threads] [-i input.gvcf.list] [-p tmp_dir]]"
    echo "Options:"
    echo "  -o STR    output file name"
    echo "  -t INT    number of threads"
    echo "  -i STR    input gvcf list"
    echo "  -p STR    temp dir (if it exists, it will be try to use the avaiable files in the dir)"
    echo "  -h        help"
    exit 1
}

while getopts "o:t:i:p:h" opt; do
    case $opt in
        o) output=$OPTARG;;
        t) threads=$OPTARG;;
        i) input=$OPTARG;;
        p) tmp_dir=$(realpath $OPTARG);;
        h) usage;;
        *) usage;;
    esac
done

if [ $# -eq 0 ] || [ -z "$output" ] || [ -z "$input" ]; then
    usage
fi

# check tools
for tool in ${GL_BIN} bcftools; do
    if ! command -v $tool &> /dev/null; then
        print_log "ERROR: $tool is not installed"
        exit 1
    fi
done

# check input files
is_exists $input


# temp dir
if [ -z "$tmp_dir" ]; then
    tmp_dir=$(pwd)/tmp_${pid}
else
    if [ ! -d $tmp_dir ]; then
        print_log "WARNING: temp dir not found: $tmp_dir"
    fi
fi

# tmp_dir=tmp_${pid}
mkdir -p $tmp_dir/vcfs
print_log "Temp dir: $(realpath $tmp_dir)"

# split gvcf per chromosome
print_log "Splitting gvcf per chromosome"
for gvcf in $(cat $input); do
    is_exists $gvcf
    if [ ! -f ${gvcf}.csi ]; then
        print_log "Indexing $gvcf"
        run_cmd "bcftools index --threads ${threads} -f $gvcf"
    fi
    run_cmd "bcftools index -s $gvcf | cut -f1 > ${tmp_dir}/chr.list.bak"

    # check if check.list between the loop are identical
    if [ -f ${tmp_dir}/chr.list ]; then
        if ! cmp -s ${tmp_dir}/chr.list ${tmp_dir}/chr.list.bak; then
            print_log "ERROR: chr.list are not identical"
            exit 1
        fi
    fi
    mv ${tmp_dir}/chr.list.bak ${tmp_dir}/chr.list

    for chr in $(cat ${tmp_dir}/chr.list); do
        if [ "$chr" == *"|"* ]; then
            continue
        fi
        if [ -f ${tmp_dir}/${gvcf/.g.vcf.gz}__${chr}.g.vcf.gz ]; then
            continue
        fi
        echo "bcftools view -r $chr $gvcf -O z -o ${tmp_dir}/${gvcf/.g.vcf.gz}__${chr}.g.vcf.gz"
    done > ${tmp_dir}/${gvcf/.g.vcf.gz}.split.sh

    if [ $(num_lines ${tmp_dir}/${gvcf/.g.vcf.gz}.split.sh) -gt 0 ]; then
        print_log "Splitting $gvcf"
        # run split gvcf
        run_cmd "parallel -j $threads < ${tmp_dir}/${gvcf/.g.vcf.gz}.split.sh"
    else
        print_log "Skipping Splitting $gvcf"
        continue
    fi
    
done

for chr in $(cat ${tmp_dir}/chr.list); do
    if [ "$chr" == *"|"* ]; then
        continue
    fi
    if [ -f ${tmp_dir}/vcfs/${chr}.vcf.gz ]; then
        continue
    fi
    ls ${tmp_dir}/*__${chr}.g.vcf.gz > ${tmp_dir}/${chr}.list
    echo "${GL_BIN} --dir ${tmp_dir}/${chr} --mem-gbytes ${threads} --config DeepVariantWGS --threads 2 --list ${tmp_dir}/${chr}.list 2> ${tmp_dir}/${chr}.log | bcftools view -O z -o ${tmp_dir}/vcfs/${chr}.vcf.gz --threads 2 > ${tmp_dir}/${chr}.log 2>&1"
done > ${tmp_dir}/glnexus.sh

if [ $(num_lines ${tmp_dir}/glnexus.sh) -gt 0 ]; then
    print_log "Running GLnexus"
    # run glnexus
    run_cmd "parallel -j $threads < ${tmp_dir}/glnexus.sh"
else
    print_log "Skipping GLnexus"
fi

# concat
print_log "Concatenating"
run_cmd "bcftools concat -o $output -O z ${tmp_dir}/vcfs/*.vcf.gz --threads $threads 2> ${tmp_dir}/concat.log"

# remove temp dir
print_log "Removing temp dir"
run_cmd "rm -rf $tmp_dir"

print_log "Done"
end_time=$(date +%s)
print_log "Elapsed time: $(seconds_to_hhmmss $((end_time - start_time)))"

#EOF