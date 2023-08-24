#!/usr/bin/env python3
import sys
import os
import rapidgzip

def count_reads(fastq_file):
    """Count the number of reads in a fastq file."""
    count = 0
    bases = 0
    line_count = 0
    with rapidgzip.open( fastq_file, parallelization = os.cpu_count() ) as file:
        file.seek( 123 )
        data = file.read( 100 )
        print( data )
            


            
        # for line in f:
        #     line_count += 1
        #     if line_count == 2:
        #         count += 1
        #         bases += len(line.strip())
        #     elif line_count == 4:
        #         line_count = 0
    return count, bases

if __name__ == '__main__':
    fastq_file = sys.argv[1]
    count, bases = count_reads(fastq_file)
    print(f'{count} reads, {bases} bases')