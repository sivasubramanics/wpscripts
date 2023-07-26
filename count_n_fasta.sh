#!/bin/bash
# Usage: count_n_fasta.sh <fasta file 1> <fasta file 2> ... <fasta file n>
# Description: Count the number of sequences in a fasta file
# Output: STDOUT (number of sequences)
# Date: 2023-07-20
# contact: siva.selvanayagam[at]wur.nl

# Check if the input file exists
c=1
for i in $@; do
    if [ $c -eq 1 ]; then
        echo -e "file_name\tn_seqs"
    fi
    if [ ! -f $i ]; then
        echo -e "$i\t0"
        continue
    fi
    c=$((c+1))
    count=$(grep -c "^>" $i)
    echo -e "$i\t$count"
done | column -t

