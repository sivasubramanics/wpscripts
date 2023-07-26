# wp-scripts
Repo containign the scripts used in the work packages. The purpose and example commandling to run such scripts are provided below.

## build_gfa.sh
This script takes a transcripts fasta file and cluster and build the gfa per clusters (considering as gene gfa).
```bash
bash build_gfa.sh transcripts.fasta output_prefix n_threads
``` 

## build_reference.sh
This script build reference resources from reference fasta and gtf/gff files. This script will merge all the contigs more than 9 chromosomes, as chromosome 0
```bash
bash build_reference.sh reference.fasta reference.gff no_chr no_threads
```

## count_reads.sh
Script to count number of reads from all the gzipped fastq file from the directory (recursively). It is recommended to install [rapidgzip](https://github.com/mxmlnkn/rapidgzip) to obtain the maximum performance
```bash
bash count_reads.sh in_dir
```

## from_yoda.py
Script to download fastq raw data from lettuceknow-release. The input meta data can be taken from the [sequencing update file](https://wageningenur4-my.sharepoint.com/:x:/r/personal/siva_selvanayagam_wur_nl/_layouts/15/Doc.aspx?sourcedoc=%7BB52F91B4-EBAD-4F31-83FE-6D368616B779%7D&file=SequencingDataMetadata_April_2023.xlsx&action=default&mobileredirect=true). Make sure you are running this in the server where already iRods is initialized (smith has it).
```bash
python from_yoda.py -m LK087.tsv -d LK087 -t download -y /nluu6p/home/research-lettuceknow-releases/1_data-releases/data-release_V1_20220921/

usage: from_yoda.py [-h] -m META_FILE -d DEST_DIR [-y YODA_DIR] -t {download,list} [-p THREADS]

Get the raw data from YODA

optional arguments:
  -h, --help            show this help message and exit
  -m META_FILE, --meta_file META_FILE
                        meta file
  -d DEST_DIR, --dest_dir DEST_DIR
                        destination directory
  -y YODA_DIR, --yoda_dir YODA_DIR
                        YODA entry directory
  -t {download,list}, --task {download,list}
                        task to perform: download, check
  -p THREADS, --threads THREADS
                        number of threads
```

## IBSpy_call_variation.sh
Script to call kmer variation from fastq file by counting kmers using [kmc](https://github.com/refresh-bio/KMC) + [kmergwas](https://github.com/voichek/kmersGWAS) and variations using [IBSpy](https://github.com/Uauy-Lab/IBSpy).
```bash
bash IBSpy_call_variation.sh -f LK001.fqlist -r genome.fa -o LK001 -p 12 -k 31 -w 50000

Usage: IBSpy_call_variation.sh [-h] [-f fq.list] [-r ref.fa] [-o output] [-p threads] [-k kmer_size] [-w window_size]
Options:
  -f STR    list of fastq files
  -r STR    reference genome
  -o STR    output file
  -p INT    number of threads
  -k INT    kmer size
  -w INT    window size
  -h        help
```

## isFile.sh
Script to return if the file exits or not
```bash
bash isFile.sh file.txt
```

## simulate_reads.sh
Script to simulate illumina paired end reads from model organism. This script downloads genome from ncbi and simulate WGS/RNA reads for a specified coverage.
```bash
bash simulate_reads.sh -t RNA -r yeast -o out -c 1 -n 2

Usage: simulate_reads.sh [-h] [-t type] [-r ref] [-o out_dir] [-p threads]
Options:
  -t STR    type of reads (WGS or RNA)
  -r STR    reference genome name (one of: 'yeast', 'fruitfly', 'human', 'mouse', 'arabidopsis', 'rice', 'zebrafish', 'maize', 'tomato', 'lettuce')
  -o STR    output directory
  -n INT    number of samples
  -c FLT    coverage
  -h        help
```
