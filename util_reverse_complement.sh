#!/bin/bash

# Function to display usage information
function usage {
    echo "Usage: $0 <input_file>"
    echo "Converts a multi-line FASTA file to its reverse complement, maintaining the same line length."
    echo "Example: $0 input.fasta"
    exit 1
}

# Function to reverse and complement a sequence
function reverse_complement {
    echo $1 | rev | tr 'ATGCatgc' 'TACGtacg'
}

# Check if user asked for help
if [[ ( $1 == "--help") ||  $1 == "-h" ]]
then 
    usage
fi

# Check if input file was provided
if [[ $# -eq 0 ]]
then
    echo "Error: No input file provided."
    usage
fi

# Initialize state and sequence accumulator
reading_sequence=false
sequence=""

# Get the length of the longest line in the fasta file
max_length=$(awk 'BEGIN{max = 0}{if($0!~/^>/){if(length$($0) > max){max = length($0);}}} END {print max}' $1)

# Reading the FASTA file
while IFS= read -r line
do
    # Check if line is an identifier
    if [[ $line =~ ^\>.* ]]
    then
        if $reading_sequence
        then
            reverse_complement $sequence | fold -w $max_length
            sequence=""
        fi
        echo $line
        reading_sequence=true
    else
        sequence+=$line
    fi
done < "$1"

# Reverse complement the last sequence
if $reading_sequence
then
    reverse_complement $sequence | fold -w $max_length
fi
