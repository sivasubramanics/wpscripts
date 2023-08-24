#!/usr/bin/env python3
# script to prepare InterPro table
# usage: python3 prepare_InterPro_table.py
# contact: siva.subramani[at]wur.nl
# date: 23-08-2023
# version: v0.1

import sys
import os

def uncomress_gz(file):
    """
    uncompress gz file
    """
    print(f"{file} processing")
    os.system(f"pigz -kd {file}")

def main():
    # wget https://ftp.ebi.ac.uk/pub/databases/interpro/current_release/interpro.xml.gz
    xml_file = "/lustre/BIF/nobackup/selva001/work/dbs/InterPro/interpro.xml.gz"
    out_file = "/lustre/BIF/nobackup/selva001/work/dbs/InterPro/interpro.tsv"

    uncomress_gz(xml_file)
    uncompressed_xml = xml_file.replace('.gz', '')

    print(f"{uncompressed_xml} processing")
    fo = open(out_file, "w")
    fo.write(f"# ipr_id\tipr_type\tipr_shortname\tipr_name\n")
    with open(uncompressed_xml, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith("<interpro id="):
                line = line.split()
                ipr_id = line[1].replace("id=", "").replace('"', '')
                ipr_shortname = line[3].replace("short_name=", "").replace('"', '')
                ipr_type = line[4].replace("type=", "").replace('"', '').replace('>', '')
                continue
            elif line.startswith("<name>"):
                ipr_name = line.replace("<name>", "").replace("</name>", "")
                continue
            elif line.startswith("</interpro>"):
                fo.write(f"{ipr_id}\t{ipr_type}\t{ipr_shortname}\t{ipr_name}\n")
                continue
            else:
                continue
    fo.close()
    print(f"{out_file} created")
    os.remove(uncompressed_xml)
    print(f"{uncompressed_xml} removed")

if __name__ == "__main__":
    main()