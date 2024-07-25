#!/bin/bash
# usage: runSTAR.sh [-h] [-f read_1.fq.gz] [-r read_2.fq.gz] [-o output_dir] [-g genome_dir] [-s genome.fasta] [-t threads]

set -e
set -u

start=$(date +%s)

# default parameters
THREADS=2
GENOME_DIR=""
OUT_DIR=""
READ1=""
READ2=""
GENOME_FASTA=""

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

usage(){
    echo "Usage: runSTAR.sh [-h] [-f read_1.fq.gz] [-r read_2.fq.gz] [-o output_dir] [-g genome_dir] [-s genome.fasta] [-t threads]"
    echo "Options:"
    echo "  -f STR    read 1 fastq file (gzipped)"
    echo "  -r STR    read 2 fastq file (gzipped)"
    echo "  -o STR    output directory"
    echo "  -g STR    STAR genome directory"
    echo "  -s STR    genome fasta file"
    echo "  -t INT    number of threads"
    echo "  -h        help"
    exit 1
}

# parse arguments
while getopts ":hf:r:o:g:t:s:" opt; do
    case $opt in
        h)
            usage
            exit 1
            ;;
        f)
            READ1=$OPTARG
            ;;
        r)
            READ2=$OPTARG
            ;;
        o)
            OUT_DIR=$OPTARG
            ;;
        g)
            GENOME_DIR=$OPTARG
            ;;
        t)
            THREADS=$OPTARG
            ;;
        s)
            GENOME_FASTA=$OPTARG
            ;;
        \?)
            echo "Invalid option: $OPTARG" 1>&2
            exit 1
            ;;
        :)
            echo "Option -$OPTARG requires an argument." 1>&2
            exit 1
            ;;
    esac
done

# check if required arguments are provided
if [ -z "$READ1" ] || [ -z "$READ2" ] || [ -z "$OUT_DIR" ] || [ -z "$GENOME_DIR" ] || [ -z "$GENOME_FASTA" ]; then
    echo "Error: missing required arguments" 1>&2
    usage
    exit 1
fi

# check if required tools are installed
for tool in STAR samtools; do
    if ! command -v $tool &> /dev/null; then
        case $tool in
            STAR) 
                echo "Error: $tool is not installed. Install STAR using the following command:" 1>&2
                echo "conda install -c bioconda star" 1>&2
                ;;
            samtools)
                echo "Error: $tool is not installed. Install samtools using the following command:" 1>&2
                echo "conda install -c bioconda samtools" 1>&2
                ;;
        esac
        exit 1
    fi
done

if [ -f "$OUT_DIR/pass2/Aligned.sortedByCoord.out.bam.bai.ok" ]; then
    print_log "STAR already completed"
    exit 0
fi

# create output directory
mkdir -p $OUT_DIR/pass1 $OUT_DIR/pass2 $OUT_DIR/genome_idx

if [ ! -f "$OUT_DIR/pass1/Aligned.out.bam.ok" ]; then
    print_log "Running STAR pass 1"
    run_cmd "STAR --runThreadN $THREADS --genomeDir $GENOME_DIR --readFilesIn $READ1 $READ2 --readFilesCommand zcat --outFileNamePrefix $OUT_DIR/pass1/ --outSAMtype BAM Unsorted > $OUT_DIR/pass1.log 2>&1"
    touch $OUT_DIR/pass1/Aligned.out.bam.ok
else
    print_log "STAR pass 1 already completed"
fi

if [ ! -f "$OUT_DIR/genome_idx/pass1.index.ok" ]; then
    print_log "Filter low coverage junctions"
    run_cmd "cat $OUT_DIR/pass1/SJ.out.tab | awk '{{ if (\$7 >= 3) print \$0}}' | sort -k1,1 -k2,2n | uniq > $OUT_DIR/pass1/junctions.txt 2> $OUT_DIR/pass1/filter_junctions.log"
    print_log "Regenerate genome index"
    run_cmd "STAR --runThreadN $THREADS --runMode genomeGenerate --genomeDir $OUT_DIR/genome_idx --genomeFastaFiles $GENOME_FASTA --sjdbFileChrStartEnd $OUT_DIR/pass1/junctions.txt --sjdbOverhang 100 > $OUT_DIR/genome_idx.log 2>&1"
    run_cmd "mv $OUT_DIR/pass1/junctions.txt $OUT_DIR/pass2/pass1_junctions.txt"
    touch $OUT_DIR/genome_idx/pass1.index.ok
else
    print_log "Genome index already generated"
fi 

if [ ! -f "$OUT_DIR/pass2/Aligned.sortedByCoord.out.bam.ok" ]; then
    print_log "Running STAR pass 2"
    run_cmd "STAR --runThreadN $THREADS --genomeDir $OUT_DIR/genome_idx --readFilesIn $READ1 $READ2 --readFilesCommand zcat --outFileNamePrefix $OUT_DIR/pass2/ --outSAMtype BAM SortedByCoordinate > $OUT_DIR/pass2.log 2>&1"
    touch $OUT_DIR/pass2/Aligned.sortedByCoord.out.bam.ok
else
    print_log "STAR pass 2 already completed"
fi

if [ ! -f "$OUT_DIR/pass2/Aligned.sortedByCoord.out.bam.bai.ok" ]; then
    print_log "Indexing BAM file"
    run_cmd "samtools index -@ $THREADS $OUT_DIR/pass2/Aligned.sortedByCoord.out.bam"
else
    print_log "BAM file already indexed"
fi

print_log "Removing intermediate files"
run_cmd "rm -rf $OUT_DIR/pass1 $OUT_DIR/genome_idx"

print_log "Done!"

end=$(date +%s)
runtime=$((end-start))
echo "Runtime: $(seconds_to_hhmmss $runtime)"


