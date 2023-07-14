#!/usr/bin/env python3
# script: from_yoda.py
# contact: siva.selvanayagam at wur.nl
# description: get the raw data from YODA
# usage: python from_yoda.py -m meta_file -d dest_dir -y yoda_dir -t task
# example: python from_yoda.py -m meta_file -d dest_dir -y yoda_dir -t download
# version: v1.0
# version date: 2023-07-13

import sys
import argparse
import os
import pandas as pd
import subprocess
import uuid

fastq_ext = ['_R1_001.fastq.gz', '_R2_001.fastq.gz', '.tar']

def is_iinit_done():
	home_dir = os.path.expanduser("~")
	irods_file = os.path.join(home_dir, ".irods", ".irodsA")
	return os.path.exists(irods_file)

def get_meta_data(meta_file):
	"""
	Get the meta data from the meta file
	:param meta_file: meta file
	:return: meta data
	"""
	meta_df = pd.read_csv(meta_file, sep="\t")
	return meta_df

def get_file(iget_path, dest_dir, threads):
	"""
	Get a file from YODA
	:param iget_path: path to file on YODA
	:param dest_dir: destination directory
	:return: None
	"""
	cmd = f"iget -N {threads} {iget_path} {dest_dir}"
	# print(f"CMD: {cmd}")
	os.system(cmd)

def list_dir(ils_path):
	"""
	List the contents of a directory on YODA
	:param ils_path: path to directory on YODA
	:return: list of files in the directory
	"""
	cmd = f"ils {ils_path}"
	filename = str(uuid.uuid4()) + ".ls"
	# print(f"CMD: {cmd} > {filename}")
	os.system(cmd + " > " + filename)
	return filename

def get_options():
	"""
	Get options from the command line
	:return: options
	"""
	parser = argparse.ArgumentParser(description='Get the raw data from YODA')
	parser.add_argument('-m', '--meta_file', type=str, help='meta file', required=True)
	parser.add_argument('-d', '--dest_dir', type=str, help='destination directory', required=True)
	parser.add_argument('-y', '--yoda_dir', type=str, default="/nluu6p/home/research-lettuceknow-releases/1_data-releases/data-release_V1_20220921/", help='YODA entry directory' , required=False)
	parser.add_argument('-t', '--task', type=str, choices = ['download', 'list'], help='task to perform: download, check', required=True)
	parser.add_argument('-p', '--threads', type=int, default=8, help='number of threads', required=False)
	args = parser.parse_args()
	return args

def check_ext(filename, ext_list):
	"""
	Check if the filename has the extension in the list
	:param filename: filename
	:param ext_list: list of extensions
	:return: True if the filename has the extension in the list, False otherwise
	"""
	for ext in ext_list:
		if ext in filename:
			return True
	return False

def main():
	
	args = get_options()
	
	if is_iinit_done():
		print(f"YODA connected")
	else:
		print(f"YODA not connected")
		sys.exit(1)
	
	# get the meta data
	meta_df = get_meta_data(args.meta_file)
	
	# get the list of accessions
	sample_ids = meta_df['SampleID_submission'].tolist()
	raw_data_dir = meta_df['Path_to_raw_data'].tolist()
	sample_seq_ids = meta_df['SampleID_Seq'].tolist()
	sequencing_run = meta_df['Sequencing_run'].tolist()
	
	# get the list of fastq files
	for idx, sample_dir in enumerate(raw_data_dir):
		dest_dir = os.path.join(args.dest_dir, sample_ids[idx], sequencing_run[idx])
		os.makedirs(dest_dir, exist_ok=True)
		list_file = list_dir(os.path.join(args.yoda_dir, sample_dir))
		fq_list = []
		with open(list_file, 'r') as f:
			for line in f:
				if sample_seq_ids[idx] in line and check_ext(line, fastq_ext):
					if args.task == 'download':
						get_file(os.path.join(args.yoda_dir, sample_dir, line.strip()), dest_dir, threads=args.threads)
					else:
						fq_list.append(line.strip())
		if args.task == 'list':
			print(f"{sample_ids[idx]}\t{len(fq_list)}\t{','.join(fq_list)}")
		os.remove(list_file)


if __name__ == "__main__":
	main()
	
