#!/usr/bin/env python3
# script to perform QC on the isoform fasta file upon on full length based on diamond search
# steps involves:
# 1. check for tools if its installed. taxonkit, diamond
# 2. download the taxonkit database if its not present
# 3. download uniref90 database if its not present
# 4. filter the uniref90 dmd database for the list of taxids
# 5. run diamond blastx
# 6. parse the diamond output
# 7. generate the summary


import argparse
import logging
import sys
import shutil
import os
import subprocess
import rapidgzip as rapidgzip

URL_UNIREF90 = "https://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref90/uniref90.fasta.gz"
URL_UNIREF90_META = "https://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref90/RELEASE.metalink"
URL_TAXONKIT = "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz"
URL_TAXONKIT_MD5 = "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz.md5"


class FASTA:
    """
    Class to store the fasta sequence
    """
    def __init__(self, name: str, sequence: str, description: str = None, readcounts: float = None):
        self.name = name
        self.sequence = sequence
        self.description = description
        self.taxid = None
        if self.description:
            if 'TaxID=' in self.description:
                self.taxid = self.description.split('TaxID=')[1].split()[0]

    def __str__(self):
        fold_seq = self.fold(60)
        return fold_seq

    def fold(self, width: int) -> str:
        out_str = f">{self.name}"
        if self.description:
            out_str += f" {self.description}"
        out_str += "\n"
        out_str += '\n'.join([self.sequence[i:i+width] for i in range(0, len(self.sequence), width)])
        return out_str

    def __len__(self):
        return len(self.sequence)


def parse_fasta(fasta_file) -> None:
    """
    Parse the fasta file
    """
    with open(fasta_file, 'r') as f:
        name = ""
        sequence = ""
        description = ""
        begun = False
        for line in f:
            if line.startswith(">"):
                if begun:
                    yield FASTA(name, sequence, description)
                line = line.strip()
                name = line[1:]
                description = line[1:]
                sequence = ""
                begun = True
            else:
                sequence += line.strip()
        yield FASTA(name, sequence, description)


def parse_fasta_gz(fasta_file, nthreads=2) -> None:
    """
    Parse the fasta file
    """

    with rapidgzip.open(fasta_file, parallelization=nthreads) as f:
        name = ""
        sequence = ""
        description = ""
        begun = False
        for line in f:
            line = line.decode('utf-8').strip()
            if line.startswith(">"):
                if begun:
                    yield FASTA(name, sequence, description)
                line = line.strip()
                name = line[1:]
                description = line[1:]
                sequence = ""
                begun = True
            else:
                sequence += line.strip()
        yield FASTA(name, sequence, description)


def run_cmd(cmd, message=None) -> list:
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


def download_file(out_file, url=None) -> None:
    """
    Download the file from the url and save it to out_file
    """
    file_name = os.path.basename(out_file)
    if file_name == 'taxdump.tar.gz':
        url = URL_TAXONKIT
        db_type = 'taxonkit'
    elif file_name == 'uniref90.fasta.gz':
        url = URL_UNIREF90
        db_type = 'uniref90'
    else:
        logging.error(f"Unknown file {file_name}")
        sys.exit(1)

    if os.path.exists(out_file):
        logging.info(f"{file_name} is already downloaded")
        if checksum_file(db_type, os.path.dirname(out_file)):
            logging.info(f"{file_name} is downloaded and the checksum is correct")
            return
        else:
            logging.error(f"{file_name} is downloaded but the checksum is incorrect. Removing file and downloading.")
            os.remove(out_file)
            download_file(out_file)
    else:
        run_cmd(f"mkdir -p {os.path.dirname(out_file)}")
    wget(url, out_file)
    if checksum_file(db_type, os.path.dirname(out_file)):
        logging.info(f"{file_name} is downloaded and the checksum is correct")
    else:
        logging.error(f"{file_name} is downloaded but the checksum is incorrect. Removing the file and retrying.")
        os.remove(out_file)
        download_file(out_file)


def check_tools() -> None:
    """
    Check for the tools if its installed: taxonkit, diamond
    """
    tools = ['taxonkit', 'diamond', 'wget', 'seqkit']
    try:
        for tool in tools:
            shutil.which(tool)
    except Exception as e:
        logging.error(f"{tool} is not installed. Please install it and add it to the PATH")
        sys.exit(1)


def wget(url, out_file) -> None:
    """
    Download the file from the url and save it to out_file
    """
    run_cmd(f"wget -q {url} -O {out_file}")


def checksum_file(dbtype, db_dir) -> bool:
    """
    Check the md5 checksum of the local file and the remote file and return True if its same else False
    """
    if dbtype == 'taxonkit':
        md5_file = os.path.join(db_dir, 'taxdump.tar.gz.md5')
        wget(URL_TAXONKIT_MD5, md5_file)
        md5_url = run_cmd(f"cat {md5_file}")[0].split()[0]
        db_local = os.path.join(db_dir, 'taxdump.tar.gz')
        md5_local = run_cmd(f"md5sum {db_local}")[0].split()[0]
        if md5_local == md5_url:
            return True
        else:
            return False
    elif dbtype == 'uniref90':
        md5_file = os.path.join(db_dir, 'RELEASE.metalink')
        wget(URL_UNIREF90_META, md5_file)
        md5_url = run_cmd(f"grep -A 4 'name=\"uniref90.fasta.gz\"' {md5_file} | grep 'type' | cut -f2 -d '>' | cut -f1 -d '<'")[0]
        db_local = os.path.join(db_dir, 'uniref90.fasta.gz')
        md5_local = run_cmd(f"md5sum {db_local}")[0].split()[0]
        if md5_local == md5_url:
            return True
        else:
            return False
    else:
        logging.error(f"Unknown dbtype {dbtype}")
        sys.exit(1)


def download_db(db_dir) -> None:
    """
    Download the taxonkit database if its not present
    """
    if not os.path.exists(db_dir):
        os.makedirs(db_dir)
    taxondir = os.path.join(db_dir, 'taxonkit')
    taxonfile = os.path.join(taxondir, 'taxdump.tar.gz')
    download_file(taxonfile)
    logging.info("Downloading taxonkit database complete.")
    # extract the taxonkit database
    run_cmd(f"tar -xzf {taxonfile} -C {taxondir}", "Extracting taxonkit database")
    logging.info("Extracting taxonkit database complete.")

    unirefdir = os.path.join(db_dir, 'uniref90')
    unireffile = os.path.join(unirefdir, 'uniref90.fasta.gz')
    download_file(unireffile)
    logging.info("Downloading uniref90 database complete.")


def get_headers(unireffile, outfile) -> int:
    """
    extract sequence ID and Taxid from the uniref90 fasta file
    """
    no_seqs = 0
    run_cmd(f"grep '^>' {unireffile} > {unireffile}.headers")
    logging.info(f"Extracting sequence ID and Taxid from the uniref90 fasta file")
    with open(f"{unireffile}.headers", 'r') as fh, open(outfile, 'w') as out:
        for line in fh:
            line = line.strip().split()
            seqid = line[0][1:]
            taxid = [x.split('=')[1] for x in line if 'TaxID' in x][0]
            if not taxid:
                logging.error(f"Taxid not found for the sequence {seqid}")
                sys.exit(1)
            out.write(f"{seqid}\t{taxid}\n")
            no_seqs += 1
    os.remove(f"{unireffile}.headers")
    return no_seqs


def get_seq_ids(unireftaxon, db_taxon, db_seqids) -> int:
    """
    Get the sequence ids for the list of taxids present in db_taxon file
    """
    no_seqs = 0
    logging.info(f"Extracting sequence IDs for the taxids")
    with open(db_taxon, 'r') as fh, open(unireftaxon, 'r') as uniref, open(db_seqids, 'w') as out:
        taxids = set([line.strip() for line in fh])
        for line in uniref:
            line = line.strip().split()
            seqid = line[0]
            taxid = line[1]
            if taxid in taxids:
                out.write(f"{seqid}\n")
                no_seqs += 1
    return no_seqs


def extract_seqs(unireffile, db_taxon, db_fasta, nthreads=2) -> int:
    """
    extract the sequences for the list of taxids
    """
    no_seqs = 0
    taxonids = []
    with open(db_taxon, 'r') as f:
        for line in f:
            line = line.strip()
            if line == '':
                continue
            if line not in taxonids:
                taxonids.append(line)
    logging.info(f"Number of taxids to be included in the DB: {len(taxonids)}")

    logging.info(f"Extracting sequences for the taxids")
    with open(db_fasta, 'w') as out:
        for fasta in parse_fasta(unireffile):
            if fasta.taxid in taxonids:
                out.write(f"{str(fasta)}\n")
                no_seqs += 1
    return no_seqs


def is_done(fname) -> None:
    """
    create a ok file with md5 checksum in it
    """
    fsize = os.path.getsize(fname)
    run_cmd(f"echo {fsize} > {fname}.ok")


def is_completed(fname) -> bool:
    """
    check if the process is complete
    """
    if os.path.exists(f"{fname}.ok"):
        fsize = os.path.getsize(fname)
        with open(f"{fname}.ok", 'r') as f:
            if fsize == float(f.readline().strip()):
                return True
    return False


def prepare_db(db_dir, dbname, nthreads, taxids=None) -> None:
    """
    Prepare the database for diamond search
    """
    taxondir = os.path.join(db_dir, 'taxonkit')
    taxonfile = os.path.join(taxondir, 'taxdump.tar.gz')
    if not os.path.exists(taxonfile):
        logging.error(f"Taxonkit database {taxonfile} is not present. Please download it.")
        logging.error("Run the command: wp1_isoform_QC.py download -d <db_dir>")
        sys.exit(1)

    # check if uniref90 fasta file is present
    unirefdir = os.path.join(db_dir, 'uniref90')
    unireffile_gz = os.path.join(unirefdir, 'uniref90.fasta.gz')
    unireffile = os.path.join(unirefdir, 'uniref90.fasta')
    if not os.path.exists(unireffile_gz):
        logging.error(f"Uniref90 database {unireffile_gz} is not present. Please download it.")
        logging.error("Run the command: wp1_isoform_QC.py download -d <db_dir>")
        sys.exit(1)

    # filter the uniref90 database for the list of taxids
    if taxids:
        os.makedirs(os.path.join(db_dir, dbname), exist_ok=True)

        # filenames
        db_fasta = os.path.join(db_dir, dbname, f"{dbname}.fasta")
        db_taxon = os.path.join(db_dir, dbname, f"{dbname}.taxids")
        db_seqids = os.path.join(db_dir, dbname, f"{dbname}.seqids")
        unireffile = os.path.join(unirefdir, 'uniref90.fasta')
        unireffile_gz = os.path.join(unirefdir, 'uniref90.fasta.gz')
        unirefseq2taxon = os.path.join(unirefdir, 'uniref90.fasta.seqids.tsv')
        if not is_completed(db_fasta):
            db_taxon = os.path.join(db_dir, dbname, f"{dbname}.taxids")
            if not is_completed(db_taxon):
                # get the children taxon ids from the taxids
                taxids = ','.join([str(taxid) for taxid in taxids])
                run_cmd(f"taxonkit list --data-dir {taxondir} --threads {nthreads} --ids {taxids} --indent '' | sort | uniq > {db_taxon}")
                is_done(db_taxon)

            if not is_completed(unireffile):
                # extract the uniref90 fasta file
                run_cmd(f"rapidgzip -kd --force -P {nthreads} {unireffile_gz}")
                is_done(unireffile)
            if not is_completed(unirefseq2taxon):
                # parse the uniref90 fasta file to get the sequence ids and its corresponding taxids
                no_seqs = get_headers(unireffile, unirefseq2taxon)
                is_done(unirefseq2taxon)
            if not is_completed(db_seqids):
                no_seqs = get_seq_ids(unirefseq2taxon, db_taxon, db_seqids)
                is_done(db_seqids)

            no_seqs = run_cmd(f"wc -l {db_seqids}")[0].split()[0]
            run_cmd(f"seqkit grep -f {db_seqids} {unireffile} > {db_fasta}")
            logging.info(f"Number of sequences to be included in the DB: {no_seqs:,}")
            is_done(db_fasta)
        # prepare the diamond database
        run_cmd(f"diamond makedb --in {db_fasta} -d {os.path.join(db_dir, dbname)}/{dbname}.dmnd -p {   nthreads}")
        logging.info("Database preparation complete.")
    else:
        run_cmd(f"diamond makedb --in {unireffile} -d {unirefdir}/uniref90")
        logging.info("Database preparation complete.")


def parse_args() -> argparse.Namespace:
    """
    Parse command line arguments
    """
    parser = argparse.ArgumentParser(
        description='This is a utility script to perform QC on the isoform fasta file upon on full length based on diamond search')
    subparsers = parser.add_subparsers(help='sub-command help', dest='command')

    # add subcommand for downloading taxonkit database and uniref90 database
    download_parser = subparsers.add_parser('download', help='download taxonkit database and uniref90 database',
                                            description='download taxonkit database and uniref90 database')
    download_parser.add_argument('-d', '--db_dir', help='database directory', required=True)

    # add subcommand for preparing the database
    prepare_parser = subparsers.add_parser('prepare', help='prepare the database for diamond search',
                                            description='prepare the database for diamond search')
    prepare_parser.add_argument('-d', '--db_dir', help='database directory', required=True)
    prepare_parser.add_argument('-n', '--name', help='database name', required=True)
    prepare_parser.add_argument('-t', '--taxids', help='taxid(s)', required=False, nargs='+', type=int)
    prepare_parser.add_argument('-p', '--threads', help='number of threads', required=False, type=int, default=2)

    # add subcommand for running diamond search and classification
    full_length_parser = subparsers.add_parser('full_length', help='run diamond search summarize full length isoforms',
                                                description='run diamond search summarize full length isoforms')
    full_length_parser.add_argument('-i', '--input', help='input fasta file', required=True)
    full_length_parser.add_argument('-o', '--output', help='output file prefix', required=True)
    full_length_parser.add_argument('-d', '--db_dir', help='diamond database directory (from prepare step)', required=True)
    full_length_parser.add_argument('-n', '--name', help='database name', required=True)
    full_length_parser.add_argument('-p', '--threads', help='number of threads', required=False, type=int, default=2)
    full_length_parser.add_argument('-f', '--force', help='force to run the command', required=False, action='store_true')

    diamond_parser = subparsers.add_parser('diamond', help='run diamond search',
                                             description='run diamond search')
    diamond_parser.add_argument('-i', '--input', help='input fasta file', required=True)
    diamond_parser.add_argument('-o', '--output', help='output file prefix', required=True)
    diamond_parser.add_argument('-d', '--db_dir', help='diamond database directory (from prepare step). basedir name will be taken from -n', required=True)
    diamond_parser.add_argument('-n', '--name', help='database name eg: "asterids" for asterids/asterids.dmnd', required=True)
    diamond_parser.add_argument('-p', '--threads', help='number of threads', required=False, type=int, default=2)
    diamond_parser.add_argument('-f', '--force', help='force to run the command', required=False, action='store_true')

    args = parser.parse_args(args=(sys.argv[1:] or ['--help']))
    return args


def full_length_summarize(input, out_prefix, db_dir, name, nthreads, is_force) -> None:
    """
    run diamond and parse the output to calculate full length summary for isoforms
    """
    diamond_out = run_diamond(db_dir, input, is_force, name, nthreads, out_prefix)

    diamond_cov = f"{out_prefix}.dmnd.cov.tsv"
    # parse the diamond output
    with open(diamond_out, 'r') as f, open(diamond_cov, 'w') as out:
        for line in f:
            line = line.strip().split()
            qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore, qcov, scov, qlen, slen = line
            if float(scov) >= 60:
                out.write(f"{qseqid}\t{qlen}\t{sseqid}\t{slen}\t{pident}\t{length}\t{qcov}\t{scov}\n")
    is_done(diamond_cov)
    logging.info(f"Full length isoform summary is written to {diamond_cov}")


def run_diamond(db_dir, input, is_force, name, nthreads, out_prefix) -> str:
    diamond_out = f"{out_prefix}.dmnd.out"
    if not is_completed(diamond_out) or is_force:
        if not os.path.exists(os.path.join(db_dir, name, f'{name}.dmnd')):
            logging.error(f"Diamond database {name} is not present. Please prepare it using 'prepare' plugin.")
            sys.exit(1)
        run_cmd(f"diamond "
                f"blastx "
                f"-d {os.path.join(db_dir, name, f'{name}.dmnd')} "
                f"-q {input} "
                f"-o {diamond_out} "
                f"-p {nthreads} "
                # f"--outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen "
                f"--sensitive --iterate -e 0.0001 --top 5  --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovhsp scovhsp qlen slen "
                f"--quiet")
        is_done(diamond_out)
    return diamond_out


def main() -> None:
    args = parse_args()
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s',
                        datefmt='%Y-%m-%d %H:%M:%S')
    check_tools()

    if args.command == 'download':
        download_db(args.db_dir)
    elif args.command == 'prepare':
        args.taxids = [int(taxid) for taxid in args.taxids]
        prepare_db(args.db_dir, args.name, args.threads, args.taxids)
    elif args.command == 'full_length':
        full_length_summarize(args.input, args.output, args.db_dir, args.name, args.threads, args.force)
    elif args.command == 'diamond':
        run_diamond(args.db_dir, args.input, args.force, args.name, args.threads, args.output)
    else:
        logging.error("Unknown command. Please check the help message.")
        sys.exit(1)


if __name__ == "__main__":
    main()
