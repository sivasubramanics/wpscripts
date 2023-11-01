#!/usr/bin/env python3
# script: prepare_pfam_table.py
# contact: siva.selvanayagam[at]wur.nl
# description: prepare pfam table


import sys
import os

def uncomress_gz(file):
    """
    uncompress gz file
    """
    print(f"{file} uncompressing")
    os.system(f"pigz -kd {file}")


def main():
    # wget https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz
    hmm_file = "/lustre/BIF/nobackup/selva001/work/dbs/Pfam/Pfam-A.hmm.gz"
    out_file = "/lustre/BIF/nobackup/selva001/work/dbs/Pfam/Pfam-A.hmm.tsv"

    uncomress_gz(hmm_file)
    uncompressed_hmm = hmm_file.replace('.gz', '')

    print(f"{uncompressed_hmm} processing")
    fo = open(out_file, "w")
    fo.write(f"# pfam_acc\tpfam_name\tpfam_desc\n")
    with open(uncompressed_hmm, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith("NAME"):
                line = line.split()
                pfam_name = line[1]
            elif line.startswith("ACC"):
                line = line.split()
                pfam_acc = line[1].split(".")[0]
            elif line.startswith("DESC"):
                line = line.split()
                pfam_desc = " ".join(line[1:])
            elif line.startswith("//"):
                fo.write(f"{pfam_acc}\t{pfam_name}\t{pfam_desc}\n")
                pfa_name = ""
                pfam_acc = ""
                pfam_desc = ""
            else:
                continue
    fo.close()
    print(f"{out_file} created")
    os.remove(uncompressed_hmm)
    print(f"{uncompressed_hmm} removed")

if __name__ == "__main__":
    main()