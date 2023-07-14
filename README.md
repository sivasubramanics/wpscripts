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
```

## IBSpy_call_variation.sh
Script to call kmer variation from fastq file by counting kmers using [kmc](https://github.com/refresh-bio/KMC) + [kmergwas](https://github.com/voichek/kmersGWAS) and variations using [IBSpy](https://github.com/Uauy-Lab/IBSpy).
```bash
bash IBSpy_call_variation.sh -f LK001.fqlist -r genome.fa -o LK001 -p 12 -k 31 -w 50000
```

## isFile.sh
Script to return if the file exits or not
```bash
bash isFile.sh file.txt
```

