#!/usr/bin/env python3

import argparse
import sys


class GTFline:
    """
    Class object to hold GTF line
    """
    def __init__(self, line):
        self.line = line.strip().split('\t')
        self.seqname = self.line[0]
        self.source = self.line[1]
        self.feature = self.line[2]
        self.start = int(self.line[3])
        self.end = int(self.line[4])
        self.score = self.line[5]
        self.strand = self.line[6]
        self.frame = self.line[7]
        self.attribute = self.line[8]
        self.attribute_dict = {}
        for attr in self.attribute.split(';'):
            if attr:
                key, value = attr.split()
                self.attribute_dict[key] = value.strip('"')
        self.length = self.end - self.start + 1

    def __str__(self):
        return '\t'.join(self.line)


def parse_gtf(gtf_file, feature):
    """
    Parse a GTF file and yield GTFline objects with the feature type specified
    :param gtf_file: GTF file
    :param feature: feature type to parse (eg: gene, transcript, exon, CDS, start_codon, stop_codon)
    :return: GTFline object
    """
    with open(gtf_file, 'r') as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            gtf = GTFline(line)
            if gtf.feature == feature:
                yield gtf


def main():
    # in_gtf = "/Users/selva001/projects/work/wp1/trinity_evigene/best_isoforms/on_expression/refmapping_SAT.lsatv11_gffcompare.annotated.gtf"
    # out_table = "/Users/selva001/projects/work/wp1/trinity_evigene/best_isoforms/on_expression/refmapping_SAT.lsatv11_gffcompare.annotated.tsv"

    parser = argparse.ArgumentParser(description="Convert cuffcompare annotated GTF to table")
    parser.add_argument("-g", "--gtf",
                        help="cuffcompare annotated GTF file (eg: refmapping_SAT.lsatv11_gffcompare.annotated.gtf)",
                        required=True, type=str)
    parser.add_argument("-o", "--out", help="Output table", required=True, type=str)
    args = parser.parse_args()

    try:
        out_fh = open(args.out, 'w')
    except IOError:
        print(f"Error: cannot open {args.out} for writing")
        sys.exit(1)

    out_fh.write("transcript_id\tcmp_ref\tclass_code\tseqname\tstart\tend\tstrand\tlength\n")
    for tr in parse_gtf(args.gtf, 'transcript'):
        transcript_id = "-"
        cmp_ref = "-"
        class_code = "-"
        if 'transcript_id' in tr.attribute_dict:
            transcript_id = tr.attribute_dict['transcript_id']
        if 'cmp_ref' in tr.attribute_dict:
            cmp_ref = tr.attribute_dict['cmp_ref']
        if 'class_code' in tr.attribute_dict:
            class_code = tr.attribute_dict['class_code']
        out_fh.write(f"{transcript_id}\t"
                     f"{cmp_ref}\t"
                     f"{class_code}\t"
                     f"{tr.seqname}\t"
                     f"{tr.start}\t"
                     f"{tr.end}\t"
                     f"{tr.strand}\t"
                     f"{tr.length}\n")


if __name__ == "__main__":
    main()