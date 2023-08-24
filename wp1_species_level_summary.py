#!/usr/bin/env python3
import sys
import os


def parse_blat(blat_psl, map_cov_cutoff, query_list=None, chr_dict=None):
    """
    parse blat psl file and return a dictionary with query_id as key and a list of target_id:target_start-target_end:map_len as value
    """
    print(f"{blat_psl} processing")
    # fo = open(f"{blat_psl}.filtered", "w")
    blat_dict = {}
    with open(blat_psl, 'r') as f:
        for iline in f:
            iline = iline.strip()
            line = iline.split()
            query_id = line[9]
            query_len = line[10]
            target_id = line[13]
            if chr_dict:
                if target_id in chr_dict:
                    target_id = chr_dict[target_id]
            target_start = line[15]
            target_end = line[16]
            map_len = line[0]
            if query_list:
                if query_id not in query_list:
                    continue
            # fo.write(f"{iline}\n")
            map_cov = float(map_len) / float(query_len)
            if map_cov >= map_cov_cutoff:
                if query_id not in blat_dict:
                    blat_dict[query_id] = []
                blat_dict[query_id].append(f"{target_id}:{target_start}-{target_end}:{int(round(map_cov, 2)*100)}")
    return blat_dict

def parse_gene_list(in_gene_list):
    """
    parse gene list file and return a list of gene ids
    """
    print(f"{in_gene_list} processing")
    gene_list = []
    with open(in_gene_list, 'r') as f:
        for line in f:
            line = line.strip().split()
            if line[0] != "gene_id":
                gene_list.append(line[0])
    return gene_list

def parse_chr_names(in_chr_names):
    """
    parse chr names file and return a dictionary with chr name as key and a list of gene ids as value
    """
    print(f"{in_chr_names} processing")
    chr_dict = {}
    with open(in_chr_names, 'r') as f:
        for line in f:
            line = line.strip().split()
            old_name = line[0]
            new_name = line[1]
            if old_name not in chr_dict:
                chr_dict[old_name] = new_name
            else:
                print(f"ERROR: {old_name} already in chr_dict")
                sys.exit(1)
    return chr_dict

def parse_pfam_hmm(in_pfam_hmm):
    print(f"{in_pfam_hmm} processing")
    pfam_hmm_dict = {}
    pfam_acc = ""
    pfam_name = ""
    pfam_desc = ""
    with open(in_pfam_hmm, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith("NAME"):
                line = line.split()
                pfam_name = line[1]
            elif line.startswith("ACC"):
                line = line.split()
                pfam_acc = line[1].split(".")[0]
                if pfam_acc not in pfam_hmm_dict:
                    pfam_hmm_dict[pfam_acc] = {}
                    pfam_hmm_dict[pfam_acc]['name'] = ""
                    pfam_hmm_dict[pfam_acc]['desc'] = ""
            elif line.startswith("DESC"):
                line = line.split()
                pfam_desc = " ".join(line[1:])
            elif line.startswith("//"):
                pfam_hmm_dict[pfam_acc]['name'] = pfam_name
                pfam_hmm_dict[pfam_acc]['desc'] = pfam_desc
                pfam_acc = ""
                pfam_name = ""
                pfam_desc = ""
            else:
                continue
    return pfam_hmm_dict

def parse_interpro(in_interpro):
    print(f"{in_interpro} processing")
    interpro_dict = {}
    with open(in_interpro, 'r') as f:
        for line in f:
            line = line.strip().split()
            ipr_id = line[0]
            ipr_type = line[1]
            ipr_shortname = line[2]
            ipr_name = " ".join(line[3:])
            if ipr_id not in interpro_dict:
                interpro_dict[ipr_id] = {}
                interpro_dict[ipr_id]['type'] = ipr_type
                interpro_dict[ipr_id]['shortname'] = ipr_shortname
                interpro_dict[ipr_id]['name'] = ipr_name
            else:
                print(f"ERROR: {ipr_id} already in interpro_dict")
                sys.exit(1)
    return interpro_dict


class Annotation:
    def __init__(self, seq_id):
        self.seq_id = seq_id
        self.length = 0
        self.pfam_acc = []
        self.pfam_name = []
        self.pfam_desc = []
        self.interpro_id = []
        self.interpro_type = []
        self.interpro_shortname = []
        self.interpro_name = []

    

def parse_annotations(in_annotations, annotations_dict, pfam_hmm_dict, interpro_dict=None):
    """
    parse annotations file and return a dictionary with gene id as key and a list of annotations as value
    """
    print(f"{in_annotations} processing")
    with open(in_annotations, 'r') as f:
        for line in f:
            line = line.strip().split('\t')
            gene_id = line[0]
            length = line[2]
            if gene_id not in annotations_dict:
                continue
            if line[3] == "Pfam":
                pfam_acc = line[4]
                if pfam_acc not in annotations_dict[gene_id].pfam_acc:
                    annotations_dict[gene_id].pfam_acc.append(pfam_acc)
                    annotations_dict[gene_id].pfam_name.append(pfam_hmm_dict[pfam_acc]['name'])
                    annotations_dict[gene_id].pfam_desc.append(pfam_hmm_dict[pfam_acc]['desc'])
            if line[11] != "-":
                interpro_id = line[11]
                if interpro_id not in annotations_dict[gene_id].interpro_id:
                    if interpro_id not in interpro_dict:
                        continue
                    annotations_dict[gene_id].interpro_id.append(interpro_id)
                    annotations_dict[gene_id].interpro_type.append(interpro_dict[interpro_id]['type'])
                    annotations_dict[gene_id].interpro_shortname.append(interpro_dict[interpro_id]['shortname'])
                    annotations_dict[gene_id].interpro_name.append(interpro_dict[interpro_id]['name'])

    return annotations_dict


def parse_fai(in_fai, seq_list=None):
    """
    parse fai file and return a dictionary with seq_id as key and length as value
    """
    print(f"{in_fai} processing")
    fo = open(f"{in_fai}.filtered", "w")
    annotations_dict = {}
    with open(in_fai, 'r') as f:
        for line in f:
            line = line.strip().split()
            seq_id = line[0]
            length = line[1]
            if seq_list:
                if seq_id not in seq_list:
                    continue
            fo.write(f"{line[0]}\t{line[1]}\n")
            if seq_id not in annotations_dict:
                annotations_dict[seq_id] = Annotation(seq_id)
            annotations_dict[seq_id].length = length
    return annotations_dict

def parse_dge(in_dge, gene_list=None):
    """
    parse dge file and return a dictionary with gene id as key and a list of dge values as value
    """
    print(f"{in_dge} processing")
    dge_dict = {}
    with open(in_dge, 'r') as f:
        for line in f:
            line = line.strip().split()
            gene_id = line[0]
            contrast = line[1]
            if contrast != "treatment":
                continue
            if gene_list:
                if gene_id not in gene_list:
                    continue
            sample = line[2][:-1]
            log2fc = round(float(line[7]),2)
            if gene_id not in dge_dict:
                dge_dict[gene_id] = {}
            if sample not in dge_dict[gene_id]:
                dge_dict[gene_id][sample] = {}
            dge_dict[gene_id][sample] = log2fc
    return dge_dict

def parse_fasta(in_fasta, seq_list=None):
    """
    parse fasta file and return a dictionary with seq_id as key and sequence as value
    """
    print(f"{in_fasta} processing")
    annotations_dict = {}
    with open(in_fasta, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                seq_id = line[1:].split()[0]
                if seq_list:
                    if seq_id not in seq_list:
                        continue
                if seq_id not in annotations_dict:
                    annotations_dict[seq_id] = Annotation(seq_id)
            else:
                if seq_id not in annotations_dict:
                    continue
                annotations_dict[seq_id].length += len(line)
    return annotations_dict

def main():
    samples = ["SAT8", "SAT12", "SAT24", "IL8", "IL12", "IL24", "SAL8", "SAL12", "SAL24"]
    chr_names = "/lustre/BIF/nobackup/selva001/work/wp1/compile_data/chr.names"
    pfam_hmm = "/lustre/BIF/nobackup/selva001/work/dbs/Pfam-A.hmm"
    interpro = "/lustre/BIF/nobackup/selva001/work/dbs/InterPro/interpro.tsv"
    
    if len(sys.argv) != 2:
        print("Usage: python3 wp1_species_level_summary.py <species>")
        print("species: IL SAL SAT all_in_one lsatv11 lsal")
        sys.exit(1)

    # specific for IL
    if sys.argv[1] == "IL":
        blat_lsal = "/lustre/BIF/nobackup/selva001/work/wp1/compile_data/kraken2/IL_trinity.okay.cont.fasta_vs_lsal_pg.psl"
        blat_lsat = "/lustre/BIF/nobackup/selva001/work/wp1/compile_data/kraken2/IL_trinity.okay.cont.fasta_vs_lsat_pg.psl"
        blat_lsatv11 = "/lustre/BIF/nobackup/selva001/work/wp1/compile_data/kraken2/IL_trinity.okay.cont.fasta_vs_lsatv11.psl"
        fasta = "/lustre/BIF/nobackup/selva001/work/wp1/compile_data/kraken2/IL_trinity.okay.cont.fasta"
        in_gene_list = "/lustre/BIF/nobackup/selva001/work/wp1/compile_data/IL.ids"
        annotations = "/lustre/BIF/nobackup/selva001/work/wp1/annotations/all_species_annotations/evigene_il.fasta.tsv"
        dge = "/lustre/BIF/nobackup/selva001/work/wp1/compile_data/IL_deg.tsv"
        out_file = "/lustre/BIF/nobackup/selva001/work/wp1/compile_data/IL.results.tsv"
    if sys.argv[1] == "SAL":
        blat_lsal = "/lustre/BIF/nobackup/selva001/work/wp1/compile_data/kraken2/SAL_trinity.okay.cont.fasta_vs_lsal_pg.psl"
        blat_lsat = "/lustre/BIF/nobackup/selva001/work/wp1/compile_data/kraken2/SAL_trinity.okay.cont.fasta_vs_lsat_pg.psl"
        blat_lsatv11 = "/lustre/BIF/nobackup/selva001/work/wp1/compile_data/kraken2/SAL_trinity.okay.cont.fasta_vs_lsatv11.psl"
        fasta = "/lustre/BIF/nobackup/selva001/work/wp1/compile_data/kraken2/SAL_trinity.okay.cont.fasta"
        in_gene_list = "/lustre/BIF/nobackup/selva001/work/wp1/compile_data/SAL.ids"
        annotations = "/lustre/BIF/nobackup/selva001/work/wp1/annotations/all_species_annotations/evigene_sal.fasta.tsv"
        dge = "/lustre/BIF/nobackup/selva001/work/wp1/compile_data/SAL_deg.tsv"
        out_file = "/lustre/BIF/nobackup/selva001/work/wp1/compile_data/SAL.results.tsv"
    if sys.argv[1] == "SAT":
        blat_lsal = "/lustre/BIF/nobackup/selva001/work/wp1/compile_data/kraken2/SAT_trinity.okay.cont.fasta_vs_lsal_pg.psl"
        blat_lsat = "/lustre/BIF/nobackup/selva001/work/wp1/compile_data/kraken2/SAT_trinity.okay.cont.fasta_vs_lsat_pg.psl"
        blat_lsatv11 = "/lustre/BIF/nobackup/selva001/work/wp1/compile_data/kraken2/SAT_trinity.okay.cont.fasta_vs_lsatv11.psl"
        fasta = "/lustre/BIF/nobackup/selva001/work/wp1/compile_data/kraken2/SAT_trinity.okay.cont.fasta"
        in_gene_list = "/lustre/BIF/nobackup/selva001/work/wp1/compile_data/SAT.ids"
        annotations = "/lustre/BIF/nobackup/selva001/work/wp1/annotations/all_species_annotations/evigene_sat.fasta.tsv"
        dge = "/lustre/BIF/nobackup/selva001/work/wp1/compile_data/SAT_deg.tsv"
        out_file = "/lustre/BIF/nobackup/selva001/work/wp1/compile_data/SAT.results.tsv"
    if sys.argv[1] == "all_in_one":
        blat_lsal = "/lustre/BIF/nobackup/selva001/work/wp1/compile_data/all_in_one_kraken2/all_in_one.cont.fasta_vs_lsal_pg.psl"
        blat_lsat = "/lustre/BIF/nobackup/selva001/work/wp1/compile_data/all_in_one_kraken2/all_in_one.cont.fasta_vs_lsat_pg.psl"
        blat_lsatv11 = "/lustre/BIF/nobackup/selva001/work/wp1/compile_data/all_in_one_kraken2/all_in_one.cont.fasta_vs_lsatv11.psl"
        fasta = "/lustre/BIF/nobackup/selva001/work/wp1/compile_data/all_in_one_kraken2/all_in_one.cont.fasta"
        in_gene_list = "/lustre/BIF/nobackup/selva001/work/wp1/compile_data/all_in_one.ids"
        annotations = "/lustre/BIF/nobackup/selva001/work/wp1/annotations/all_species_annotations/evigene_all.fasta.tsv"
        dge = "/lustre/BIF/nobackup/selva001/work/wp1/compile_data/all_in_one_deg.tsv"
        out_file = "/lustre/BIF/nobackup/selva001/work/wp1/compile_data/all_in_one.results.tsv"
    if sys.argv[1] == "lsatv11":
        blat_lsal = "/lustre/BIF/nobackup/selva001/genomes/LSat_v11_ncbi/ncbi_dataset/data/GCF_002870075.4/rna.fna_vs_lsal_pg.psl"
        blat_lsat = "/lustre/BIF/nobackup/selva001/genomes/LSat_v11_ncbi/ncbi_dataset/data/GCF_002870075.4/rna.fna_vs_lsat_pg.psl"
        blat_lsatv11 = "/lustre/BIF/nobackup/selva001/genomes/LSat_v11_ncbi/ncbi_dataset/data/GCF_002870075.4/rna.fna_vs_lsatv11.psl"
        fasta = "/lustre/BIF/nobackup/selva001/genomes/LSat_v11_ncbi/ncbi_dataset/data/GCF_002870075.4/rna.fna"
        in_gene_list = "/lustre/BIF/nobackup/selva001/work/wp1/compile_data/lsat.ids"
        annotations = "/lustre/BIF/nobackup/selva001/work/wp1/annotations/all_species_annotations/lsat_mod.fasta.tsv"
        dge = "/lustre/BIF/nobackup/selva001/work/wp1/compile_data/lsat_deg.tsv"
        out_file = "/lustre/BIF/nobackup/selva001/work/wp1/compile_data/lsat.results.tsv"
    if sys.argv[1] == "lsal":
        blat_lsal = "/lustre/BIF/nobackup/selva001/genomes/LSal/resources/Lsal_rnd2.all.maker.transcripts.renamed.fasta_vs_lsal_pg.psl"
        blat_lsat = "/lustre/BIF/nobackup/selva001/genomes/LSal/resources/Lsal_rnd2.all.maker.transcripts.renamed.fasta_vs_lsat_pg.psl"
        blat_lsatv11 = "/lustre/BIF/nobackup/selva001/genomes/LSal/resources/Lsal_rnd2.all.maker.transcripts.renamed.fasta_vs_lsatv11.psl"
        fasta = "/lustre/BIF/nobackup/selva001/genomes/LSal/resources/Lsal_rnd2.all.maker.transcripts.renamed.fasta"
        in_gene_list = "/lustre/BIF/nobackup/selva001/work/wp1/compile_data/lsal.ids"
        annotations = "/lustre/BIF/nobackup/selva001/work/wp1/annotations/all_species_annotations/lsal.fasta.tsv"
        dge = "/lustre/BIF/nobackup/selva001/work/wp1/compile_data/lsal_deg.tsv"
        out_file = "/lustre/BIF/nobackup/selva001/work/wp1/compile_data/lsal.results.tsv"


    gene_list = parse_gene_list(in_gene_list)
    annotations_dict = parse_fasta(fasta, gene_list)
    chr_dict = parse_chr_names(chr_names)
    blat_lsal_dict = parse_blat(blat_lsal, 0.8, gene_list, chr_dict)
    blat_lsat_dict = parse_blat(blat_lsat, 0.8, gene_list, chr_dict)
    blat_lsatv11_dict = parse_blat(blat_lsatv11, 0.8, gene_list, chr_dict)
    pfam_hmm_dict = parse_pfam_hmm(pfam_hmm)
    interpro_dict = parse_interpro(interpro)
    annotations_dict = parse_annotations(annotations, annotations_dict, pfam_hmm_dict, interpro_dict)
    dge_dict = parse_dge(dge, gene_list)

    print(f"number of genes before decont: {len(gene_list)}")
    gene_list = list(set(gene_list) & set(annotations_dict.keys()))
    print(f"number of genes after decont: {len(gene_list)}")

    fo = open(out_file, "w")
    fo.write(f"gene_id\tlength\tpfam_name\tinterpro_shortname\tlsal_pg\tlsat_pg\tlsatv11")
    for sample in samples:
        fo.write(f"\t{sample}")
    fo.write(f"\n")
    for gene in gene_list:
        fo.write(f"{gene}")
        fo.write(f"\t{annotations_dict[gene].length}")
        # COLUMN 3: PFAM NAME
        if annotations_dict[gene].pfam_name:
            fo.write(f"\t{join_elements(annotations_dict[gene].pfam_name)}")
        else:
            fo.write(f"\t-")

        # COLUMN 4: INTERPRO SHORTNAME
        if annotations_dict[gene].interpro_shortname:
            fo.write(f"\t{join_elements(annotations_dict[gene].interpro_shortname)}")
        else:
            fo.write(f"\t-")

        # COLUMN 5: LSAL_PG
        if gene in blat_lsal_dict:
            fo.write(f"\t{join_elements(blat_lsal_dict[gene])}")
        else:
            fo.write(f"\t-")

        # COLUMN 6: LSAT_PG
        if gene in blat_lsat_dict:
            fo.write(f"\t{join_elements(blat_lsat_dict[gene])}")
        else:
            fo.write(f"\t-")

        # COLUMN 7-: LSATV11
        if gene in blat_lsatv11_dict:
            fo.write(f"\t{join_elements(blat_lsatv11_dict[gene])}")
        else:
            fo.write(f"\t-")

        # COLUMN 8: DGE
        if gene in dge_dict:
            for sample in samples:
                if sample in dge_dict[gene]:
                    fo.write(f"\t{dge_dict[gene][sample]}")
                else:
                    fo.write(f"\t-")
        else:
            for sample in samples:
                fo.write(f"\t-")
        fo.write(f"\n")
    fo.close()

def join_elements(elements):
    """
    Join elements of a list with a comma and a space.
    """
    # print(elements)
    return ", ".join(elements)

if __name__ == "__main__":
    main()
    
                






