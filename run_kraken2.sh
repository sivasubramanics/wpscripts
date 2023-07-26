#!/bin/bash
# Usage: run_kraken2.sh <input directory> <output directory> <kraken2 database> <threads>
# Description: Run kraken2 on all the fasta files in the input directory
# Output: Kraken2 output files in the output directory
# Date: 2023-07-20
# contact: siva.selvanayagam[at]wur.nl

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
    echo "Usage: $0 [-h] [-i in1.fasta,in2.fasta,...inN.fasta] [-o output_dir] [-d kraken2_db] [-t threads]"
    echo "Options:"
    echo "  -i STR    input fasta files"
    echo "  -o STR    output directory"
    echo "  -d STR    kraken2 database"
    echo "  -t INT    number of threads"
    echo "  -h        help"
    exit 1
}

# get options
while getopts "i:o:d:t:h" opt; do
    case $opt in
        i) in_fasta=$OPTARG;;
        o) out_dir=$OPTARG;;
        d) DB=$OPTARG;;
        t) THREADS=$OPTARG;;
        h) usage;;
        *) usage;;
    esac
done

# check if required tools are installed
for tool in kraken2 ktImportTaxonomy extract_kraken_reads.py; do
    if ! command -v $tool &> /dev/null; then
        case $tool in
            kraken2) 
                echo "Error: $tool is not installed. Install kraken2 using the following command:"
                echo "conda install -c bioconda kraken2"
                ;;
            ktImportTaxonomy)
                echo "Error: $tool is not installed. Install ktImportTaxonomy using the following command:"
                echo "conda install -c bioconda krona"
                ;;
            extract_kraken_reads.py)
                echo "Error: $tool is not installed. Install extract_kraken_reads.py using the following command:"
                echo "conda install -c bioconda krakentools"
                ;;
        esac
        exit 1
    fi
done

# check if required arguments are provided
if [ $# -eq 0 ] || [ -z "$in_fasta" ] || [ -z "$out_dir" ] || [ -z "$DB" ] || [ -z "$THREADS" ]; then
    usage
fi
{
    ORGDB=$DB
    TAXONS="2 6656 554915 2698737"
}

# create output directory
print_log "Creating output directory"
run_cmd "mkdir -p $out_dir"

# check if input fasta is one file or multiple files
if [[ "$in_fasta" =~ "," ]]; then
    n_inputs=$(echo $in_fasta | tr ',' '\n' | wc -l)
    in_fasta=$(echo $in_fasta | tr ',' ' ')
else
    n_inputs=1
fi

# check if kraken2 database is in /dev/shm
if [ $n_inputs -gt 1 ]; then
    size_shm=$(df /dev/shm | tail -n 1 | awk '{print $4/1048576}' | bc -l)
    print_log "Size of /dev/shm: $size_shm GB"
    size_db=$(du -s $DB | awk '{print $1/1048576}' | bc -l)
    print_log "Size of $DB: $size_db GB"
    if [ $(echo "$size_shm / $size_db" | bc) -ge 2 ]; then
        print_log "Copying $DB to /dev/shm"
        cp -r -n $DB /dev/shm/
        ORGDB=$DB
        DB=/dev/shm/$(basename $DB)
    fi
fi

# run kraken2
for in_fa in $in_fasta;do
    out_prefix=$(basename $in_fa)
    out_prefix=${out_prefix%.*}
    if [ ! -f $in_fa ]; then
        echo "Error: $in_fa does not exist"
        exit 1
    fi
    print_log "Running kraken2 on $in_fa"
    run_cmd "kraken2 --db $DB --threads $THREADS --output $out_dir/$out_prefix.kraken2.out --report $out_dir/$out_prefix.kraken2.report.txt $in_fa &> $out_dir/$out_prefix.kraken2.log"

    print_log "Extract columns 2 and 3 from kraken2.out"
    run_cmd "cut -f 2,3 $out_dir/$out_prefix.kraken2.out > $out_dir/$out_prefix.kraken2.out.krona"

    print_log "Running ktImportTaxonomy on $out_dir/$out_prefix.kraken2.out.krona"
    run_cmd "ktImportTaxonomy $out_dir/$out_prefix.kraken2.out.krona -o $out_dir/$out_prefix.kraken2.out.krona.html &> $out_dir/$out_prefix.ktImportTaxonomy.log"

    print_log "Possible contaminants in $in_fa"
    run_cmd "awk '\$1>1 && \$4==\"S\"' $out_dir/$out_prefix.kraken2.report.txt > $out_dir/$out_prefix.kraken2.contaminants"

    print_log "Running extract_kraken_reads.py on $out_dir/$out_prefix.kraken2.out"
    run_cmd "extract_kraken_reads.py -k $out_dir/$out_prefix.kraken2.out -r $out_dir/$out_prefix.kraken2.report.txt -s $in_fa --include-children -o $out_dir/$out_prefix.cont.fasta --exclude -t $TAXONS &> $out_dir/$out_prefix.extract_kraken_reads.log"
done

# remove kraken2 database from /dev/shm
if [ "$ORGDB" != "$DB" ]; then
    print_log "Removing $DB from /dev/shm"
    run_cmd "rm -rf $DB"
fi

# end clock
end=$(date +%s)
runtime=$((end-start))

print_log "Finished in $(seconds_to_hhmmss $runtime)"

# EOF