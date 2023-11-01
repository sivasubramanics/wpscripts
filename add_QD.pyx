# cython: language_level=3
import gzip
import os
import sys

cdef bint check_if_gz(str invcf):
    if invcf.endswith('.gz'):
        return True
    else:
        return False

cpdef add_QD(str vcf, str output):
    cdef int n_samples
    cdef double dp, qual, qd
    cdef str line, info, i
    cdef list line_list, info_list
    if check_if_gz(vcf):
        vcf_fh = gzip.open(vcf, 'rt')
    else:
        vcf_fh = open(vcf, 'r')
    if check_if_gz(output):
        output_fh = gzip.open(output, 'wt')
    else:
        output_fh = open(output, 'w')
    
    for line in vcf_fh:
        if line.startswith('##'):
            output_fh.write(line)
        elif line.startswith('#'):
            line_list = line.split('\t')
            n_samples = len(line_list) - 9
            output_fh.write(f"##INFO=<ID=QD,Number=1,Type=Float,Description=\"QUAL/DP\">\n")
            output_fh.write(f"##vcf_mod_command=vcf.py {' '.join(sys.argv[1:])}\n")
            output_fh.write(line)
        else:
            line = line.rstrip()
            line_list = line.split('\t')
            info = line_list[7]
            info_list = info.split(';')
            dp = 0.0
            for i in info_list:
                if i.startswith('QD='):
                    print('QD already present in VCF file')
                    os.remove(output)
                    sys.exit(1)
                if i.startswith('DP='):
                    dp = float(i.split('=')[1])
                    dp = dp / n_samples
            qual = float(line_list[5])
            qd = round(qual / dp, 2)
            info_list.append(f'QD={qd}')
            info = ';'.join(info_list)
            line_list[7] = info
            line = '\t'.join(line_list)
            print(f"{qual}:{dp}:{qd}")
            output_fh.write(line + '\n')
