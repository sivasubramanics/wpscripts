#!/usr/bin/env python3

import argparse
import gzip
import sys
import os
from Bio import bgzf


def main():
    parser = argparse.ArgumentParser(description='script to process VCF files')
    subparsers = parser.add_subparsers(dest='command', help='commands')

    add_QD_parser = subparsers.add_parser('add_QD', help='add QD to VCF file')
    add_QD_parser.add_argument('-v', '--vcf', help='input VCF file', required=True)
    add_QD_parser.add_argument('-o', '--output', help='output VCF file', required=True)
    add_QD_parser.add_argument('-t', '--threads', help='number of threads to use', default=2, type=int)

    args = parser.parse_args(args=(sys.argv[1:] or ['--help']))

    if args.command == 'add_QD':
        add_QD(args.vcf, args.output, args.threads)

def check_if_gz(invcf):
    if invcf.endswith('.gz'):
        return True
    else:
        return False

def add_QD(vcf, output, threads=2):
    if check_if_gz(vcf):
        vcf_fh = bgzf.BgzfReader(vcf, 'r')
    else:
        vcf_fh = open(vcf, 'r')
    if check_if_gz(output):
        # output_fh = gzip.open(output, 'wt')
        output_fh = open(output.replace('.gz', ''), 'w')
    else:
        output_fh = open(output, 'w')
    
    dp_flag = False
    qd_flag = False
    for line in vcf_fh:
        if line.startswith('##'):
            output_fh.write(line)
            if line.startswith('##INFO=<ID=DP'):
                dp_flag = True
            if line.startswith('##INFO=<ID=QD'):
                qd_flag = True
        elif line.startswith('#'):
            n_samples = len(line.split('\t')) - 9
            if not dp_flag:
                output_fh.write(f"##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">\n")
            if not qd_flag:
                output_fh.write(f"##INFO=<ID=QD,Number=1,Type=Float,Description=\"QUAL/mean_DP\">\n")
                output_fh.write(f"##vcf_mod_command=vcf.py {' '.join(sys.argv[1:])}\n")
            output_fh.write(line)
        else:
            line = line.rstrip()
            line = line.split('\t')
            info = line[7]
            info = info.split(';')
            if qd_flag:
                line = '\t'.join(line)
                output_fh.write(line + '\n')
                continue
            elif dp_flag and not qd_flag:
                for i in info:
                    if i.startswith('QD='):
                        print('QD already present in VCF file')
                        os.remove(output)
                        sys.exit(1)
                    if i.startswith('DP='):
                        dp = i.split('=')[1]
                        dp = int(dp)
                        mean_dp = round(dp/n_samples, 2)
            else:
                dp = 0
                try :
                    dp_idx = line[8].split(':').index('DP')
                except ValueError:
                    print('DP not present in VCF file')
                    os.remove(output.replace('.gz', ''))
                    sys.exit(1)
                for i in line[9:]:
                    dp += int(i.split(':')[dp_idx])
                mean_dp = round(dp/n_samples, 2)
                info.append(f'DP={dp}')
            qual = float(line[5])
            qd = round(qual / mean_dp, 2)
            info.append(f'QD={qd}')
            info = ';'.join(info)
            line[7] = info
            line = '\t'.join(line)
            output_fh.write(line + '\n')
    if check_if_gz(output):
        output_fh.close()
        os.system(f'bgzip --threads {threads} -f {output.replace(".gz", "")}')


if __name__ == '__main__':
    main()

