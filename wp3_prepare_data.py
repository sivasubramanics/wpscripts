#!/usr/bin/env python3
# script to prepare data for transcriptome assembly pipeline

import os
import argparse
import logging
import subprocess
import sys


def run_cmd(cmd, message=None):
    """
    Run the command and return the stdout and stderr as string list as [stdout, stderr]
    """
    try:
        if message:
            logging.info(message)
        logging.info(f"CMD: {cmd}")
        process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = process.communicate()
        stdout = stdout.decode('utf-8').strip()
        stderr = stderr.decode('utf-8').strip()
        run_code = process.returncode
        if run_code != 0:
            logging.error(f"Error running the command: {cmd}")
            logging.error(f"Error: {stderr}")
            sys.exit(1)
        return stdout, stderr
    except Exception as e:
        logging.error(f"Error running the command: {cmd}")
        logging.error(f"Error: {e}")
        sys.exit(1)

def untar(in_dir):
    """
    # extract all the tar files if it exists in any nexted directory. Make sure the tar extraction is done in the same directory where the tar file exists.
    """
    extracted_files = []
    for root, dirs, files in os.walk(in_dir):
        for file in files:
            if file.endswith('.tar'):
                tar_file = os.path.join(root, file)
                stdout, stderr = run_cmd(f'tar -xvf {tar_file} -C {os.path.dirname(tar_file)}', f"Extracting {tar_file}")
                files = stdout.split('\n')
                for f in files:
                    if f:
                        # append absolute path of extracted file
                        extracted_files.append(os.path.abspath(os.path.join(os.path.dirname(tar_file), f)))
    return extracted_files

def get_fastqs(in_dir):
    """
    Get fastq filenames from a sample directory. search recursively for fastq files, yeilding pairs of fastq files.
    yield sample_name, run_name, rep_name, fwd, rev
    Ideally the path for each fastq will be {sample_name}/{run_name}/{rep_name}_R1.fastq.gz
    """
    for root, dirs, files in os.walk(in_dir):
        for file in files:
            if file.endswith('_R1_001.fastq.gz'):
                fwd = os.path.join(root, file)
                rev = fwd.replace('_R1_001.fastq.gz', '_R2_001.fastq.gz')
                if os.path.exists(rev):
                    # ideally the path for each fastq will be {sample_name}/{run_name}/{rep_name}_R1.fastq.gz
                    acc_name = os.path.basename(os.path.dirname(os.path.dirname(os.path.dirname(fwd))))
                    sample_name = os.path.basename(os.path.dirname(os.path.dirname(fwd)))
                    run_name = os.path.basename(os.path.dirname(fwd))
                    rep_name = os.path.basename(fwd).replace('_R1_001.fastq.gz', '')
                    yield [acc_name, sample_name, run_name, rep_name, fwd, rev]


def main():
    parser = argparse.ArgumentParser(description='Prepare data for transcriptome assembly pipeline')
    parser.add_argument('-i', '--input', help='Input directory containing raw data', required=True)
    parser.add_argument('-o', '--output', help='Output directory to store processed data', required=True)
    args = parser.parse_args()

    # logging base config
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s',
                        datefmt='%Y-%m-%d %H:%M:%S')

    # create output directory
    os.makedirs(args.output, exist_ok=True)

    # extract all the tar files if it exists in any nexted directory
    extracted_file = untar(args.input)
    if extracted_file:
        with open(f"{args.output}/extracted_files.txt", 'w') as out:
            for f in extracted_file:
                out.write(f"{f}\n")

    with open(os.path.join(args.output, 'metadata.tsv'), 'w') as out:
        out.write('#accession\tsample\trep\tfq1\tfq2\n')
        for acc,sample,run,rep,fwd,rev in get_fastqs(args.input):
            if not acc:
                logging.error(f"Accession name not found for {fwd}")
                sys.exit(1)
            if not sample:
                logging.error(f"Sample name not found for {fwd}")
                sys.exit(1)
            if not run:
                logging.error(f"Run name not found for {fwd}")
                sys.exit(1)
            if not rep:
                logging.error(f"Replicate name not found for {fwd}")
                sys.exit(1)
            os.makedirs(os.path.join(args.output, acc, sample, run), exist_ok=True)
            new_fwd = os.path.join(args.output, acc, sample, run, os.path.basename(fwd))
            new_rev = os.path.join(args.output, acc, sample, run, os.path.basename(rev))
            abs_fwd = os.path.abspath(fwd)
            abs_rev = os.path.abspath(rev)
            os.symlink(abs_fwd, new_fwd)
            os.symlink(abs_rev, new_rev)
            out.write(f'{acc}\t{sample}\t{rep}\t{os.path.abspath(new_fwd)}\t{os.path.abspath(new_rev)}\n')



if __name__ == '__main__':
    main()
