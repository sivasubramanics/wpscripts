#!/bin/bash
# script: convert_format.sh
# description: convert format from different formats
# contact: siva.selvanayagam[at]wur.nl
# date: 2023-07-24
# version: 0.0.1

set -u
set -e

threads=1

print_log(){
    echo "[`date +'%Y-%m-%d %H:%M:%S'`] $1"
}

run_cmd(){
    echo "[`date +'%Y-%m-%d %H:%M:%S'`] CMD: $1"
    eval $1
}

check_tool(){
    if ! command -v $1 &> /dev/null; then
        print_log "Error: $1 is not installed"
        exit 1
    fi
}

usage(){
    echo "Usage: $0 [-h] [-i input] [-o output] [-e in_format] [-f out_format] [-t threads]"
    echo "Options:"
    echo "  -i STR    input file"
    echo "  -o STR    output file"
    echo "  -e STR    input file format"
    echo "  -f STR    output file format"
    echo "  -t INT    number of threads"
    echo "  -h        help"
    exit 1
}

while getopts "i:o:e:f:t:h" opt; do
    case $opt in
        i) input=$OPTARG;;
        o) output=$OPTARG;;
        e) in_format=$OPTARG;;
        f) out_format=$OPTARG;;
        t) threads=$OPTARG;;
        h) usage;;
        *) usage;;
    esac
done

if [ $# -eq 0 ] || [ -z "$input" ] || [ -z "$output" ] || [ -z "$in_format" ] || [ -z "$out_format" ]; then
    usage
fi

case $in_format in
    vcf)
        case $out_format in
            vcf.gz)
                check_tool bcftools
                run_cmd "bcftools view --threads $threads -Oz -o $output $input > convert_format.log 2>&1"
                ;;
            plink)
                check_tool plink
                run_cmd "plink --vcf $input --make-bed --out $output > convert_format.log 2>&1"
                ;;
            hmp)
                check_tool run_pipeline.pl
                run_cmd "run_pipeline.pl -fork1 -vcf $input -sortPositions -export $output -exportType Hapmap -runfork1 > convert_format.log 2>&1"
                ;;
            *)
                print_log "Error: $out_format is not supported"
                exit 1
                ;;
        esac
        ;;
    vcf.gz)
        case $out_format in
            vcf)
                check_tool bcftools
                run_cmd "bcftools view --threads $threads -Ov -o $output $input > convert_format.log 2>&1"
                ;;
            plink)
                check_tool plink
                run_cmd "plink --vcf $input --make-bed --out $output > convert_format.log 2>&1"
                ;;
            hmp)
                check_tool run_pipeline.pl
                run_cmd "run_pipeline.pl -fork1 -vcf $input -sortPositions -export $output -exportType Hapmap -runfork1 > convert_format.log 2>&1"
                ;;
            *)
                print_log "Error: $out_format is not supported"
                exit 1
                ;;
        esac
        ;;
    bam)
        case $out_format in
            sam)
                check_tool samtools
                run_cmd "samtools view -@ $threads -o $output $input > convert_format.log 2>&1"
                ;;
            bed)
                check_tool bedtools
                run_cmd "bedtools bamtobed -i $input > $output 2> convert_format.log"
                ;;
            *)
                print_log "Error: $out_format is not supported"
                exit 1
                ;;
        esac
        ;;
    sam)
        case $out_format in
            bam)
                check_tool samtools
                run_cmd "samtools view -@ $threads -b -o $output $input > convert_format.log 2>&1"
                ;;
            bed)
                check_tool bedtools
                run_cmd "bedtools bamtobed -i $input > $output 2> convert_format.log"
                ;;
            *)
                print_log "Error: $out_format is not supported"
                exit 1
                ;;
        esac
        ;;
    *)
        print_log "Error: $in_format is not supported"
        exit 1
        ;;
esac