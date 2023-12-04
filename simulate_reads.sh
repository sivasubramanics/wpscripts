#!/bin/bash
# script: simulate_reads.sh
# description: simulate reads from reference genome
# contact: siva.selvanayagam[at]wur.nl
# date: 2023-07-19
# version: 0.1

set -u
set -e

# start clock
start=$(date +%s)

# defaults for formatting
YELLOW=$'\e[1;33m'
GREEN=$'\e[0;32m'
RED=$'\e[0;31m'
NC=$'\e[0m'
ULINE=$'\e[4m'

# default arguments
n_samples=1
n_reps=3
read_length=100
insert_len=700
insert_std=200
mapping=false
threads=1

# function to write fasta_length.sh
write_fasta_length_script(){
    echo -e "#!/bin/bash"
    echo -e "set -u"
    echo -e "set -e"
    echo -e "c=1"
    echo -e "for i in \$@; do"
    echo -e "    if [ ! -f \$i ]; then"
    echo -e "       echo \"Error: \$i does not exist\""
    echo -e "       exit 1"
    echo -e "   fi"
    echo -e "   if [ \$c -eq 1 ]\; then"
    echo -e "       echo -e \"file_name\tlength\""
    echo -e "   fi"
    echo -e "   c=\$((c+1))"
    echo -e "   awk '/^>/{if (seqlen){print substr(seqid, 2), seqlen;};seqid = \$1;seqlen = 0; next; }{seqlen += length(\$0)}END {print substr(seqid, 2), seqlen}' \$i > \$i.length"
    echo -e "   awk -v filename=\$i '{s+=\$2} END {print filename, s}' \$i.length"
    echo -e "done | column -t"
}


# function to print log
print_log(){
    echo "[`date +'%Y-%m-%d %H:%M:%S'`] LOG: ${YELLOW}${ULINE}$1${NC}"
}


# function to run command
run_cmd(){
    echo "[`date +'%Y-%m-%d %H:%M:%S'`] CMD: $1"
    eval $1
}


# function to get coverage
get_cov(){
    local mean=$1
    local std_dev=${2:-0.3}
    local min=$(echo "$mean - ($mean * $std_dev)" | bc)
    local max=$(echo "$mean + ($mean * $std_dev)" | bc)
    local range=$(echo "$max - $min + 1" | bc)
    local rand=$(awk 'BEGIN{srand(); print rand()}')
    echo $(echo "$min + $rand * ($max - $min)" | bc)
}


# function to get length of fasta file
get_length() {
    local file="$1"
    local length=0
    while read -r line; do
        if [[ $line == ">"* ]]; then
            echo "$line"
        else
            len=${#line}
            length=$((length + len))
            echo "$line" | tr -d '\n'
            echo " $len"
            fi
    done < "$file"
    echo "$length"
}


# function to convert seconds to human readable time
seconds_to_hhmmss() {
  local seconds=$1
  local hours=$(printf "%02d" $((seconds / 3600)))
  local minutes=$(printf "%02d" $(((seconds / 60) % 60)))
  local seconds=$(printf "%02d" $((seconds % 60)))
  echo "$hours:$minutes:$seconds"
}


# function to print options
print_options(){
    echo -e "--------------"
    echo -e ${ULINE}"INPUT OPTIONS:"${NC}
    echo -e "type\t$type"
    echo -e "ref_name\t$ref_name ($ref)"
    echo -e "out_dir\t$out_dir"
    echo -e "n_samples\t$n_samples"
    echo -e "coverage\t$coverage"
    echo -e "mapping\t$mapping"
    echo -e "threads\t$threads"
    echo -e "--------------"
}


# usage message
usage(){
    echo "Usage: $0 [-h] [-t type] [-r ref] [-o out_dir] [-p threads] [-m]"
    echo "Options:"
    echo "  -t STR    type of reads (WGS or RNA)"
    echo "  -r STR    reference genome name (one of: 'yeast', 'fruitfly', 'human', 'mouse', 'arabidopsis', 'rice', 'zebrafish', 'maize', 'tomato', 'lettuce')"
    echo "  -o STR    output directory"
    echo "  -n INT    number of samples"
    echo "  -c FLT    coverage"
    echo "  -p INT    number of threads"
    echo "  -m        mapping"
    echo "  -h        help"
    echo "-----------------------------------------------------------------------------------------"
    echo "Example     \`bash simulate_reads.sh -t WGS -r yeast -o sim_reads_out -n 2 -c 2 -p 8 -m\`"
    echo "RUN TIME    ~1 min"
    echo "-----------------------------------------------------------------------------------------"
    exit 1
}


# parse arguments
while getopts "t:r:o:n:c:p:mh" opt; do
    case $opt in
        t) type=$OPTARG;;
        r) ref_name=$OPTARG;;
        o) out_dir=$OPTARG;;
        n) n_samples=$OPTARG;;
        c) coverage=$OPTARG;;
        p) threads=$OPTARG;;
        m) mapping=true;;
        h) usage;;
        *) usage;;
    esac
done


# check arguments
if [ $# -eq 0 ] || [ -z "$type" ] || [ -z "$ref_name" ] || [ -z "$out_dir" ] || [ -z "$coverage" ]; then
    usage
fi


# check tools
for tool in datasets wgsim pigz fasta_length.sh; do
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
            wgsim)
                echo "Error: wgsim is not installed. Please install wgsim using the following command:"
                echo "conda install -c bioconda wgsim"
                ;;
            pigz)
                echo "Error: pigz is not installed. Please install pigz using the following command:"
                echo "conda install -c conda-forge pigz"
                ;;
            fasta_length.sh)
                echo "Error: fasta_length.sh is not found in \$PATH"
                echo "I am creating a file called fasta_length.sh in the current directory. Please add this file to your PATH [OR] run the following commands:"
                write_fasta_length_script > fasta_length.sh
                echo ""
                echo "mkdir -p $(realpath $(echo $HOME))/bin"
                echo "chmod +x fasta_length.sh"
                echo "mv fasta_length.sh $(realpath $(echo $HOME))/bin/"
                echo "export PATH=$(realpath $(echo $HOME))/bin:\$PATH"
                echo ""
                ;;
        esac
        exit 1
    fi
done

if [ "$mapping" == "true" ]; then
    for tool in bwa-mem2 samtools; do
        if ! command -v $tool &> /dev/null; then
            case $tool in
                bwa)
                    echo "Error: bwa is not installed. Please install bwa using the following command:"
                    echo "conda install -c bioconda bwa-mem2"
                    ;;
                samtools)
                    echo "Error: samtools is not installed. Please install samtools using the following command:"
                    echo "conda install -c bioconda samtools"
                    ;;
            esac
            exit 1
        fi
    done
fi


# check reference genome
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


# -------------------- RUN BEGINS HERE --------------------
# print options
print_options | column -t -s $'\t'

# create output directory
mkdir -p $out_dir
mkdir -p $out_dir/download
mkdir -p $out_dir/logs
mkdir -p $out_dir/genome
mkdir -p $out_dir/reads
if [ "$mapping" == "true" ]; then
    mkdir -p $out_dir/bams
fi
ref_genome_fa="$out_dir/genome/$ref_name.genome.fa"
ref_rna_fa="$out_dir/genome/$ref_name.rna.fa"
ref_gff="$out_dir/genome/$ref_name.genome.gff"


# download reference genome
if [ ! -f $ref_genome_fa ] || [ ! -f $ref_rna_fa ] || [ ! -f $ref_gff ]; then
    print_log "Downloading reference genome $ref_name"
    run_cmd "datasets download genome accession GCF_000146045.2 --include gff3,rna,genome --filename $out_dir/download/ncbi_dataset.zip --no-progressbar > $out_dir/logs/download.log 2>&1"
    print_log "Unzipping reference genome"
    run_cmd "unzip -o -q $out_dir/download/ncbi_dataset.zip -d $out_dir/download > $out_dir/logs/unzip.log 2>&1"
    print_log "Moving reference genome to $out_dir/genome"
    genome_file=$(find $out_dir/download -name *_genomic.fna)
    run_cmd "mv $genome_file $ref_genome_fa"
    run_cmd "mv $out_dir/download/ncbi_dataset/data/GCF_000146045.2/genomic.gff $ref_gff"
    run_cmd "mv $out_dir/download/ncbi_dataset/data/GCF_000146045.2/rna.fna $ref_rna_fa"
fi


# get reference file name, index and length
if [ "$type" == "WGS" ]; then
    ref_fa=$ref_genome_fa
fi
if [ "$type" == "RNA" ]; then
    ref_fa=$ref_rna_fa
fi
length=$(fasta_length.sh $ref_fa | tail -n1 | awk '{print $NF}')
if [ "$mapping" == "true" ]; then
    if [ ! -f $ref_fa.0123 ]; then
        print_log "Indexing reference genome"
        run_cmd "bwa-mem2 index $ref_fa > $out_dir/logs/index.log 2>&1"
    fi
fi


# simulate reads
for i in $(seq 1 $n_samples); do
    if [ "$type" == "WGS" ]; then
        print_log "Simulating reads for sample $i/$n_samples"
        prefix=$ref_name"_"$type"_"sample_$i
        cov=$(get_cov $coverage 0.1)
        n_reads=$(echo "$cov * $length / 100 / 2" | bc)
        run_cmd "wgsim -d $insert_len -s $insert_std -1 $read_length -2 $read_length -N $n_reads $ref_fa $out_dir/reads/$prefix.1.fq $out_dir/reads/$prefix.2.fq > $out_dir/reads/$prefix.vars 2> $out_dir/logs/wgsim_$prefix.log"
        run_cmd "pigz -f $out_dir/reads/$prefix.1.fq $out_dir/reads/$prefix.2.fq"
        if [ "$mapping" == "true" ]; then
            rm -f $out_dir/logs/mapping_$prefix.log
            print_log "Mapping reads for sample $i/$n_samples"
            run_cmd "bwa-mem2 mem -t $threads -R '@RG\tID:$prefix\tSM:$prefix' $ref_fa $out_dir/reads/$prefix.1.fq.gz $out_dir/reads/$prefix.2.fq.gz 2>> $out_dir/logs/mapping_$prefix.log | samtools sort -@ $((threads-1)) -o $out_dir/bams/$prefix.bam - >> $out_dir/logs/mapping_$prefix.log 2>&1"
            print_log "Indexing bam file for sample $i/$n_samples"
            run_cmd "samtools index -@ $threads $out_dir/bams/$prefix.bam"
        fi
    fi
    if [ "$type" == "RNA" ]; then
        cov=$(get_cov $coverage 0.5)
        for j in $(seq 1 $n_reps); do
            print_log "Simulating reads for sample $i/$n_samples - rep $j/$n_reps"
            prefix=$ref_name"_"$type"_"sample_$i"_"rep_$j
            cov=$(get_cov $cov 0.001)
            n_reads=$(echo "$cov * $length / 100 / 2" | bc)
            run_cmd "wgsim -d $insert_len -s $insert_std -1 $read_length -2 $read_length -N $n_reads $ref_fa $out_dir/reads/$prefix.1.fq $out_dir/reads/$prefix.2.fq > $out_dir/reads/$prefix.vars 2> $out_dir/logs/wgsim_$prefix.log"
            run_cmd "pigz -f $out_dir/reads/$prefix.1.fq $out_dir/reads/$prefix.2.fq"
            if [ "$mapping" == "true" ]; then
                rm -f $out_dir/logs/mapping_$prefix.log
                print_log "Mapping reads for sample $i/$n_samples - rep $j/$n_reps"
                run_cmd "bwa-mem2 mem -t $threads -R '@RG\tID:$prefix\tSM:$prefix' $ref_fa $out_dir/reads/$prefix.1.fq.gz $out_dir/reads/$prefix.2.fq.gz 2>> $out_dir/logs/mapping_$prefix.log | samtools sort -@ $((threads-1)) -o $out_dir/bams/$prefix.bam - >> $out_dir/logs/mapping_$prefix.log 2>&1"
                print_log "Indexing bam file for sample $i/$n_samples - rep $j/$n_reps"
                run_cmd "samtools index -@ $threads $out_dir/bams/$prefix.bam"
            fi
        done
    fi
done


# remove intermediate files
run_cmd "rm -rf ref $out_dir/download"
# print tree view of output directory
run_cmd "tree --noreport -h $out_dir/reads $out_dir/bams"
# --------------------- RUN ENDS HERE ---------------------

# end clock
end=$(date +%s)
runtime=$((end-start))

# print runtime
print_log "Finished in $(seconds_to_hhmmss $runtime)"
