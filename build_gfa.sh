#!/bin/bash


set -u
set -e

start=$(date +%s)

print_log(){
    echo "[`date +'%Y-%m-%d %H:%M:%S'`] $1"
}

run_cmd(){
    echo "[`date +'%Y-%m-%d %H:%M:%S'`] CMD: $1"
    eval $1
}

num_seqs(){
    grep -c "^>" $1
}

if [ $# -ne 3 ]; then
    echo "Usage: $0 <transcripts.fasta> <out_dir> <threads>"
    exit 1
fi

# check if the list of tools are available
for tool in samtools mmseqs clustalo faSomeRecords PanPA parallel; do
    if ! command -v $tool &> /dev/null; then
        echo "\"$tool\" could not be found, please install and rerun the script"
        exit
    fi
done

# input transcript fasta file
transcripts=$1
out_dir=$2
threads=$3

# check if the output directory exists, if not create it. If exists, exit
if [ -d $out_dir ]; then
    echo "$out_dir already exists"
    exit
else
    mkdir -p $out_dir
    mkdir -p $out_dir/logs
    mkdir -p $out_dir/genes
    mkdir -p $out_dir/msa
    mkdir -p $out_dir/graphs
    mkdir -p $out_dir/clusters
    mkdir -p $out_dir/sh
fi

# generate samtools index
print_log "Generating samtools index"
samtools faidx $transcripts > $out_dir/logs/samtools-faidx.log 2>&1

# run clustering using mmseqs
print_log "Running clustering using mmseqs"
mmseqs easy-linclust $transcripts $out_dir/clusters/clust_out $out_dir/clusters/tmp --threads $threads > $out_dir/logs/mmseqs.log 2>&1

# generate clusters files
print_log "Generating clusters files"
awk -v out=$out_dir/clusters/clust '{if($1!=prev){i++;print $2 > out"_"i".ids"}else{print $2 >> out"_"i".ids"};prev=$1}' $out_dir/clusters/clust_out_cluster.tsv

# generate clusters fasta files
print_log "Splitting the clusters"
for file in $out_dir/clusters/clust_*.ids; do
    echo "faSomeRecords $transcripts $file $file.fasta" >> $out_dir/sh/splitfasta.sh
done
parallel -j $threads < $out_dir/sh/splitfasta.sh > $out_dir/logs/splitfasta.log 2>&1

# move the cluster_fasta to genes directory
mv $out_dir/clusters/clust_*.ids.fasta $out_dir/genes/

# generate msa for each cluster
print_log "Generating msa for each cluster"
for file in $out_dir/genes/clust_*.ids.fasta; do
    if [ $(num_seqs $file) -lt 2 ]; then
        cp $file $file.msa.fasta
        continue
    fi
    echo "clustalo -i $file -o $file.msa.fasta --force > $file.msa.log 2>&1" >> $out_dir/sh/msa.sh
done
parallel -j $threads < $out_dir/sh/msa.sh > $out_dir/logs/msa.log 2>&1

# move the cluster_fasta.msa to msa directory
mv $out_dir/genes/clust_*.ids.fasta.msa.fasta $out_dir/msa/
mv $out_dir/genes/clust_*.ids.fasta.msa.log $out_dir/logs/

# generate graph for each cluster
print_log "Generating graph for each cluster"
PanPA build_gfa -d $out_dir/msa -o $out_dir/graphs -c $threads > $out_dir/logs/PanPA.log 2>&1

# index the graph
print_log "Indexing the graph"
PanPA build_index -d $out_dir/msa --seeding_alg wk_min -k 5 -w 3 --seed_limit 0 -o $out_dir/index_k5_w3_no_limit.index > $out_dir/logs/PanPA-build_index.log 2>&1



# calculate the time taken for the complete run
end=$(date +%s)
runtime=$((end-start))

print_log "Runtime: $runtime seconds"
