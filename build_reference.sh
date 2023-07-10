#!/bin/bash

set -x
set -e
set -u

if [ $# -ne 4 ]; then
    echo "Usage: $0 <reference fasta> <reference gff> <number of chromosomes> <threads>"
    exit 1
fi

# # reference genome fasta and gff files
# fasta=/lustre/BIF/nobackup/selva001/work/references/from_ncbi/Lactuca_sativa_v11/ncbi_dataset/data/GCF_002870075.4/GCF_002870075.4_Lsat_Salinas_v11_genomic.fna
# gff=/lustre/BIF/nobackup/selva001/work/references/from_yoda/merged-annotation-Lactuca-sativa-V11-RefSeq-Maker_V1_20230202/GCF_002870075.4_Lsat_Salinas_v11_Fused.gff
# threads=20
# n_chr=9

# reference genome fasta and gff files
fasta=$1
gff=$2
n_chr=$3
threads=$4

# generate fai for the downloaded reference fasta
samtools faidx $fasta

# generate chromosome map to covnert chromosome names to numbers
awk -v n=$n_chr '{i++;if(i>n){print $1"\t"0}else{print $1"\t"i}}' $fasta.fai > chr.map

# convert the chromosome names to numbers on fasta file and gff file
merge_contigs.py --out genome \
    --fasta $fasta \
    --gff $gff \
    --gap 200 \
    --chrmap chr.map > merge_contigs.log 2>&1

# convert gff to gtf
agat_convert_sp_gff2gtf.pl --gff genome.gff -o genome.gtf > gff2gtf.log 2>&1

# replace the chromosome names with numbers
sed 's/^SEQ/0/' -i genome.gtf 

# extract protein sequences
agat_sp_extract_sequences.pl -g genome.gtf -f genome.fa -t cds -p -o proteins.faa > genome2proteins.log 2>&1

# extract transcript sequences
agat_sp_extract_sequences.pl -g genome.gtf -f genome.fa -t exon --merge -o transcripts.fna > genome2transcripts.log 2>&1

# extract cds sequences
agat_sp_extract_sequences.pl -g genome.gtf -f genome.fa -t cds -o cds.fna > genome2cds.log 2>&1

# find identical transcripts
pblat -threads=$threads transcripts.fna transcripts.fna transcripts_selfmap.psl > pblat.log 2>&1
awk '{if($10!=$14 && $1==$11 && $1==$15){print}}' transcripts_selfmap.psl > transcripts_selfmap.identical.psl

# generate samtools index
samtools faidx genome.fa > samtools-faidx.log 2>&1

# generate minimap2 index
minimap2 -t $threads -d genome.mmi genome.fa > minimap2.mmi.log 2>&1

# generate miniprot index
miniprot -t $threads -d genome.mpi genome.fa > miniprot.log 2>&1

# generate bowtie index
bowtie-build --threads $threads genome.fa genome > bowtie-build.log 2>&1

# generate bowtie2 index
bowtie2-build --threads $threads genome.fa genome > bowtie2-build.log 2>&1

# makeblastdb
makeblastdb -in proteins.faa -dbtype prot -out proteins > makeblastdb.log 2>&1

# generate bwa index
bwa-mem2 index -p genome genome.fa > bwa-mem2.log 2>&1

# generate hisat2 index
# hisat2-build --threads $threads genome.fa genome > hisat2-build.log 2>&1

# generate salmon index
# salmon index -t transcripts.fna -i transcripts > salmon-index.log 2>&1

# generate kallisto index
# kallisto index -i transcripts.kallisto.idx transcripts.fna > kallisto-index.log 2>&1

# generate star index
# STAR --runThreadN $threads --runMode genomeGenerate --genomeDir star --genomeFastaFiles genome.fa --sjdbGTFfile genome.gtf > star-index.log 2>&1


# move the log files to logs directory
mkdir logs
mv *.log logs

