#!/usr/bin/env python3

"""
This script can be used to convert DArTseq SNP_mapping_2.csv file to HMP file.
"""

import argparse
from collections import defaultdict

nucl_dict = {'AA': 'A', 'AC': 'M', 'AG': 'R', 'AT': 'W', 'CA': 'M', 'GA': 'R', 'TA': 'W',
             'CC': 'C', 'CG': 'S', 'CT': 'Y', 'GC': 'S', 'TC': 'Y',
             'GG': 'G', 'GT': 'K', 'TG': 'K',
             'TT': 'T', 'NN': 'N'}


class SNP:
    def __init__(self, header, data, chr_col, pos_col, snp_col, strand_col, sample_col, pic_col, chrstring='Chr'):
        self.chr = data[chr_col].replace(chrstring, '')
        self.pos = data[pos_col]
        self.snp = data[snp_col]
        if data[strand_col] == 'Minus':
            self.strand = '-'
        else:
            self.strand = '+'
        self.ref_allele, self.alt_allele = self.get_alleles()
        self.gt_num = data[sample_col:]
        self.samples = header[sample_col:]
        self.gt_nucl = self.gt2alleles()
        self.alleles = self.get_major_minor_alleles()
        self.id = f'S{self.chr}_{self.pos}'
        self.clone_id = data[1]
        self.pic = data[pic_col]

    def get_major_minor_alleles(self):
        """
        from the genotype alleles, calculate major and minor alleles and return them as list
        """
        bases = []
        alleles = []
        count = {}
        for nucl in self.gt_nucl:
            if nucl == 'NN':
                continue
            bases.extend([nucl[0], nucl[1]])
        # count the number of each base and sort them based on number of occurrences
        for base in set(bases):
            count[base] = bases.count(base)
        count = sorted(count.items(), key=lambda x: x[1], reverse=True)
        for i in count:
            alleles.append(i[0])
        return alleles

    def get_alleles(self):
        """
        Get reference and alternate alleles
        """
        alleles = self.snp.split(':')
        alleles = alleles[1].split('>')
        return alleles[0], alleles[1]

    def gt2alleles(self):
        """
        Get all alleles for the genotype numbers
        """
        alleles = []
        for i in self.gt_num:
            if i == '0':
                alleles.append(f'{self.ref_allele}{self.ref_allele}')
            elif i == '1':
                alleles.append(f'{self.ref_allele}{self.alt_allele}')
            elif i == '2':
                alleles.append(f'{self.alt_allele}{self.alt_allele}')
            else:
                alleles.append('NN')
        return alleles


def parse_dart_snp_mapping_2(dart_snp_mapping_2, reference, chrstring='Chr'):
    """
    Parse DArTseq SNP_mapping_2.csv file
    """
    with open(dart_snp_mapping_2, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('*'):
                line = line.split(',')
                sample_col = 0
                for i in line:
                    if i == '*':
                        sample_col += 1
                    else:
                        break
                continue
            if line.startswith('AlleleID'):
                header = line.split(',')
                chr_col = header.index(f'Chrom_{reference}')
                pos_col = header.index(f'ChromPosSnp_{reference}')
                snp_col = header.index('SNP')
                strand_col = header.index(f'Strand_{reference}')
                pic_col = header.index('AvgPIC')
                continue
            data = line.split(',')
            snp = SNP(header, data, chr_col, pos_col, snp_col, strand_col, sample_col, pic_col, chrstring)
            yield snp


def compare_snp(snp_one, snp_two):
    """
    comparing two SNPs and see if there is any mismatch in the genotype alleles for any sample
    """
    for i in range(len(snp_one.gt_nucl)):
        if snp_one.gt_nucl[i] != snp_two.gt_nucl[i]:
            if snp_one.gt_nucl[i] == 'NN' or snp_two.gt_nucl[i] == 'NN':
                continue
            # print(f'Mismatch found for SNP: {snp_one.id} in sample: {snp_one.samples[i]}. '
            #       f'Alleles: {snp_one.samples[i]}:{snp_one.gt_nucl[i]} - {snp_two.samples[i]}:{snp_two.gt_nucl[i]}')



def main():
    parser = argparse.ArgumentParser(description='Convert DArTseq SNP_mapping_2.csv file to HMP file')
    parser.add_argument('-i', '--input', help='Input DArTseq SNP_mapping_2.csv file', required=True)
    parser.add_argument('-r', '--reference', help='Reference genome name (you may find it in the file. eg: '
                                                  'Chrom_Rice_RGAP_v7: Rice_RGAP_v7)', required=True)
    parser.add_argument('-o', '--output', help='Output HMP file', required=False)
    parser.add_argument('-d', '--two_letter', help='Two letter code for the alleles', action='store_true')
    parser.add_argument('-s', '--chrstring', help='String to be used for chromosome name', required=False)
    args = parser.parse_args()

    if not args.output:
        args.output = f'{args.input.rstrip(".csv")}.hmp.txt'

    if args.chrstring:
        chrstring = args.chrstring

    genotype_data = defaultdict(SNP)
    for snp in parse_dart_snp_mapping_2(args.input, args.reference, chrstring):
        if snp.pos == '0':
            # remove SNPs with position 0
            continue
        if snp.id in genotype_data:
            if genotype_data[snp.id].pic < snp.pic:
                genotype_data[snp.id] = snp
        else:
            genotype_data[snp.id] = snp

    with open(args.output, 'w') as f:
        nSNPs = 0
        for snp in genotype_data.values():
            nSNPs += 1
            if snp.chr == '':
                nSNPs -= 1
                continue
            if nSNPs == 1:
                f.write(f'rs#\talleles\tchrom\tpos\tstrand\tassembly#\tcenter\tprotLSID\tassayLSID\tpanelLSID\t'
                        f'QCcode\t{"\t".join(snp.samples)}\n')
            alleles = []
            for i in snp.gt_nucl:
                if args.two_letter:
                    alleles.append(i)
                else:
                    alleles.append(nucl_dict[i])
            f.write(f'{snp.id}_{snp.clone_id}\t{"/".join(snp.alleles)}\t{snp.chr}\t{snp.pos}\t{snp.strand}\tNA\tNA\tNA\tNA\tNA\tNA\t'
                    f'{"\t".join(map(str, snp.gt_nucl))}\n')


if __name__ == '__main__':
    main()
