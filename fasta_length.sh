#!/bin/bash
# script: fasta_length.sh
# description: get length of sequence in fasta file
# contact: siva.selvanayagam[at]wur.nl
# date: 2023-07-21
# version: 0.1

set -u
set -e

c=1
for i in $@; do
    if [ ! -f $i ]; then
        echo "Error: $i does not exist"
        exit 1
    fi
    if [ $c -eq 1 ]; then
        echo -e "file_name\tlength"
    fi
    c=$((c+1))
    awk '/^>/{if (seqlen){print substr(seqid, 2), seqlen;};seqid = $1;seqlen = 0; next; }{seqlen += length($0)}END {print substr(seqid, 2), seqlen}' $i > $i.length
    awk -v filename=$i 'BEGIN{OFMT="%.0f";}{s+=$2} END {print filename, s}' $i.length 
done | column -t