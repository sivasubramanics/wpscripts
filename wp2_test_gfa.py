#!/usr/bin/env python3

import sys
import os
import argparse
import time
import subprocess
import pandas as pd
import networkx as nx
from datetime import datetime
from collections import defaultdict

bufsize = 0
# start clock
start_time = time.time()
psl_head = "matches\tmismatches\trep.matches\tNs\tQgapcount\tQgapbases\tTgapcount\tTgapbases\tstrand\tQname\tQsize\tQstart\tQend\tTname\tTsize\tTstart\tTend\tblockcount\tblocksizes\tqStarts\ttStarts"
nthreads = 2


class GAF(object):
    """
    GAF object (referenced from GraphAligner output)
    """
    def __init__(self, line):
        line = line.strip().split('\t')
        self.read_name = line[0]
        self.read_length = int(line[1])
        self.read_start = int(line[2])
        self.read_end = int(line[3])
        self.strand = line[4]
        self.path_name = line[5]
        self.path_length = int(line[6])
        self.path_start = int(line[7])
        self.path_end = int(line[8])
        self.mapping_quality = int(line[9])
        self.mismatch = int(line[10])
        self.aln_length = int(line[11])
        self.tags = {}
        for i in range(12, len(line)):
            tag = line[i].split(':')
            if tag[0] not in self.tags:
                self.tags[tag[0]] = tag[2]
            else:
                print(f"Duplicate tag {tag[0]} found in GAF file")
                sys.exit(1)
    
    def nodes(self):
        """
        Get the nodes from the path
        """
        nodes = []
        directions = []
        i = 0
        for char in self.path_name:
            if char in ['<', '>']:
                directions.append(char)
                nodes.append('')
                i += 1
            else:
                nodes[i-1] += char
        return directions, nodes


class FASTA(object):
    def __init__(self, name, sequence, description=None):
        self.name = name
        self.sequence = sequence.upper()
        self.description = description
    
    def __str__(self):
        return f">{self.name}\n{self.sequence}"
    
    def __len__(self):
        return len(self.sequence)
    
    def rev_complement(self):
        return FASTA(self.name, reverse_complement(self.sequence))

    def write_seq(self, handle):
        handle.write(f">{self.name}\n{self.sequence}\n")
    
    def write_aln(self, handle, strand='+'):
        if len(self.name) > 20:
            name = self.name[:20]
        else:
            name = self.name + " "*(20-len(self.name)) 
        handle.write(f"{name}\t{strand}\t{self.sequence}\n")

class PSL(object):
    def __init__(self, line):
        self.matches = int(line[0])
        self.misMatches = int(line[1])
        self.repMatches = int(line[2])
        self.nCount = int(line[3])
        self.qNumInsert = int(line[4])
        self.qBaseInsert = int(line[5])
        self.tNumInsert = int(line[6])
        self.tBaseInsert = int(line[7])
        self.strand = line[8]
        self.qName = line[9]
        self.qSize = int(line[10])
        self.qStart = int(line[11])
        self.qEnd = int(line[12])
        self.tName = line[13]
        self.tSize = int(line[14])
        self.tStart = int(line[15])
        self.tEnd = int(line[16])
        self.blockCount = int(line[17])
        self.blockSizes = [int(x) for x in line[18].split(",") if x]
        self.qStarts = [int(x) for x in line[19].split(",") if x]
        self.tStarts = [int(x) for x in line[20].split(",") if x]
        self.qCoverage = round(sum(self.blockSizes) / self.qSize, 2)
        self.tCoverage = round(sum(self.blockSizes) / self.tSize, 2)
        self.inline = line

    def __str__(self):
        return ("\t".join(self.inline))
    
    def write_psl(self, handle):
        inline = "\t".join(self.inline)
        handle.write(f"{inline}\n")

class Segment:
    """
    Segment object (you can call it as Node as well)
    """
    def __init__(self, name, sequence):
        self.name = name
        self.sequence = sequence
        self.length = len(sequence)
    
    def __str__(self):
        return f"S\t{self.name}\t{self.sequence}\tLN:i:{self.length}"
    
    def __repr__(self):
        return str(self)
    
class Link:
    """
    Link object (you can call it as Edge as well)
    """
    def __init__(self, from_segment, from_strand, to_segment, to_strand, overlap):
        self.from_segment = from_segment
        self.to_segment = to_segment
        self.from_strand = from_strand
        self.to_strand = to_strand
        self.overlap = overlap

    def __str__(self):
        return f"L\t{self.from_segment}\t{self.from_strand}\t{self.to_segment}\t{self.to_strand}\t{self.overlap}"
    
    def __repr__(self):
        return str(self)
    
class Path:
    """
    Path object
    """
    def __init__(self, name, segments, directions, overlaps):
        self.name = name
        self.segments = segments
        self.directions = directions
        self.overlaps = overlaps
    
    def __str__(self):
        return f"P\t{self.name}\t{','.join([s+d for s, d in zip(self.segments, self.directions)])}\t{','.join(self.overlaps)}"
    
    def __repr__(self):
        return str(self)
    
class GFA:
    """
    GFA object
    May cause large memory usage if the gfa file is large
    """
    def __init__(self, segments, links, paths):
        self.segments = segments
        self.links = links
        self.paths = paths

    def __str__(self):
        return "\n".join([str(segment) for segment in self.segments] + [str(link) for link in self.links] + [str(path) for path in self.paths])
    
    def __repr__(self):
        return str(self)
    
    def summary(self):
        return (f"Number of segments: {len(self.segments)}\n"
                f"Number of links: {len(self.links)}\n"
                f"Number of paths: {len(self.paths)}")
    
    def extract_path(self, path_name, output_file):
        print(f"Extracting path: {path_name}")
        path = [p for p in self.paths if p.name == path_name][0]
        print(f"Path length: {path}")
        sequence = ''.join([s.sequence for s in self.segments if s.name in path.segments])
        with open(output_file, 'w') as file:
            file.write(f'>{path_name}\n')
            file.write(fold(sequence, 60))

    def sub_graph(self, node_name):
        paths = [p for p in self.paths if node_name in p.segments]
        segment_names = [name for path in paths for name in path.segments]
        segments = [s for s in self.segments if s.name in segment_names]
        links = [l for l in self.links if l.from_segment in segment_names or l.to_segment in segment_names]
        return GFA(segments, links, paths)
    
    def write_to_file(self, filepath):
        with open(filepath, 'w') as file:
            file.write(str(self))

class GFAParser:
    """
    GFA parser
    """
    def __init__(self, filepath):
        self.filepath = filepath

    def parse(self):
        segments = []
        links = []
        paths = []
        print(f"Parsing GFA file: {self.filepath}")
        with open(self.filepath, 'r') as file:
            for line in file:
                fields = line.strip().split("\t")
                if fields[0] == 'S':
                    segments.append(Segment(fields[1], fields[2]))
                elif fields[0] == 'L':
                    links.append(Link(fields[1], fields[2], fields[3], fields[4], fields[5]))
                elif fields[0] == 'P':
                    segs = fields[2].split(',')
                    seg_names = [s[:-1] for s in segs]
                    directions = [s[-1] for s in segs]
                    overlaps = fields[3].split(',')
                    paths.append(Path(fields[1], seg_names, directions, overlaps))
        return GFA(segments, links, paths)

class Blocks(object):
    """
    Blocks from PSL
    """
    def __init__(self, name):
        self.name = name
        self.block_seqs = []
        self.breaks = []
        self.path = []
    
    def add_breaks(self, breaks):
        self.breaks.extend(breaks)
    
    def update_block_seqs(self, sequence):
        self.breaks = list(sorted(set(self.breaks)))
        self.block_seqs = []
        if self.breaks[-1] != len(sequence):
            self.breaks.append(len(sequence))
        if self.breaks[0] != 0:
            self.breaks.insert(0, 0)
        for i in range(len(self.breaks)-1):
            seqs = sequence[self.breaks[i]:self.breaks[i+1]]
            # if len(seqs) > 256: # split into 256bp chunks since some of the gfa parsing tools doesn't accept the nodes longer than 256bp
            #     self.block_seqs.extend(split_segments(seqs))
            # else:
            #     self.block_seqs.append(seqs)
            self.block_seqs.append(seqs)

    def get_links(self, handle):
        for i in range(len(self.path)-1):
            handle.write(f"L\t{self.path[i]}\t+\t{self.path[i+1]}\t+\t0M\n")

    def write_path(self, handle):
        handle.write(f"P\t{self.name}\t{','.join([str(x)+'+' for x in self.path])}\t{','.join(['0M' for x in self.path])}\n")
    

def reverse_complement(sequence):
    """
    Get the reverse complement of a sequence
    """
    sequence = sequence.upper()
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return "".join([complement[base] for base in sequence[::-1]])

def parse_fasta(fasta_file):
    """
    Parse a fasta file
    """
    with open(fasta_file, 'r') as f:
        name = ""
        sequence = ""
        begun = False
        for line in f:
            if line.startswith(">"):
                line = line.strip()
                line = line.split(" ")
                if begun:
                    yield FASTA(name, sequence, description)
                name = line[0].replace('>', '')
                if len(line) > 1:
                    description = " ".join(line[1:])
                else:
                    description = None
                sequence = ""
                begun = True
            else:
                sequence += line.strip()
        yield FASTA(name, sequence, description)

def parse_psl(psl_file):
    """
    Parse a psl file
    """
    with open(psl_file, 'r') as f:
        for line in f:
            line = line.strip().split("\t")
            yield PSL(line)

def get_seq(sequence, blocks_sizes, block_starts):
    """
    Get the sequence from a block
    """
    sequence_len = len(sequence)
    last_block_end = 0
    seqs = []
    fseq = ""
    for nblock in range(len(blocks_sizes)):
        block_size = blocks_sizes[nblock]
        block_start = block_starts[nblock]
        block_end = block_start + block_size
        if block_start > last_block_end:
            seqs.append(sequence[last_block_end:block_start])
        if block_size == 0:
            continue
        seqs.append(sequence[block_start:block_end])
        if nblock == len(blocks_sizes) - 1 and block_end < sequence_len:
            seqs.append(sequence[block_end:])
        last_block_end = block_end
    return seqs
        
def split_segments(node_seq, chunk_size=256):
    """
    Split a sequence into 256bp chunks, since construction of vg index is limited to 256bp segments
    """
    chunks = []
    for i in range(0, len(node_seq), chunk_size):
        chunks.append(node_seq[i:i+chunk_size])
    return chunks

def get_breaks(psl_record):
    """
    Get the breaks from a block
    """
    qbreaks = [0, psl_record.qSize] + [item for sublist in [[psl_record.qStarts[i], psl_record.qStarts[i] + psl_record.blockSizes[i]] for i in range(psl_record.blockCount)] for item in sublist]
    
    if psl_record.strand == '+':
        tbreaks = [0, psl_record.tSize] + [item for sublist in [[psl_record.tStarts[i], psl_record.tStarts[i] + psl_record.blockSizes[i]] for i in range(psl_record.blockCount)] for item in sublist]
    else:
        tbreaks = [0, psl_record.tSize] + [item for sublist in [[psl_record.tSize - psl_record.tStarts[i], psl_record.tSize - psl_record.tStarts[i] - psl_record.blockSizes[i]] for i in range(psl_record.blockCount)] for item in sublist]
    
    return qbreaks, tbreaks

def psl_to_gfa(args):
    """
    Convert a psl file to a gfa file
    """
    
    if args.gfa.endswith('.gfa'):
        args.out = args.gfa.replace('.gfa', '')
    else:
        print(f"Output file must be a gfa file. '.gfa' extension is required")
        sys.exit(1)

    if args.outpsl:
        flt_psl = open(args.outpsl, 'w')

    block_dict = defaultdict()
    combinations = defaultdict(dict)
    transcript_direction = {}
    edges = open(f"{args.out}.edges.tsv", 'w')
    edges.write(f"query"
                 f"\ttarget"
                 f"\tqcoverage"
                 f"\ttcoverage"
                 f"\tnum_blocks"
                 f"\tblock_size"
                 f"\t{psl_head}\n")
    
    print_log(f"Parsing PSL file")
    for record in parse_psl(args.psl):
        if record.qName == record.tName:
            continue
        if record.qCoverage < args.coverage or record.tCoverage < args.coverage:
            continue
        if record.blockCount > args.blocks:
            continue
        if args.outpsl:
            record.write_psl(flt_psl)
        # add breaks to block_dict
        if record.qName not in block_dict:
            block_dict[record.qName] = Blocks(record.qName)
        if record.tName not in block_dict:
            block_dict[record.tName] = Blocks(record.tName)
        qbreaks, tbreaks = get_breaks(record)
        block_dict[record.qName].add_breaks(qbreaks)
        block_dict[record.tName].add_breaks(tbreaks)

        # add edges between transcripts
        if record.qName not in combinations:
            combinations[record.qName][record.tName] = record
        elif record.tName not in combinations[record.qName]:
            combinations[record.qName][record.tName] = record
        edges.write(f"{record.qName}"
                    f"\t{record.tName}"
                    f"\t{record.qCoverage}"
                    f"\t{record.tCoverage}"
                    f"\t{record.blockCount}"
                    f"\t{record.blockSizes}"
                    f"\t{record}\n")
        transcript_direction = update_strand_info(transcript_direction, record.qName, record.tName, record.strand)
    edges.close()
    
    print_log(f"Parsing Segments")
    nodes = []
    seq_dict = {}
    for transcript in parse_fasta(args.fasta):
        seq_dict[transcript.name] = transcript
        if transcript.name not in block_dict:
            block_dict[transcript.name] = Blocks(transcript.name)
            block_dict[transcript.name].add_breaks([0, len(transcript.sequence)])
        block_dict[transcript.name].update_block_seqs(transcript.sequence)
        nodes.extend(block_dict[transcript.name].block_seqs)
    nodes = index_dict(set(nodes))

    print(f"---------- psl to gfa ----------\n"
          f"Number of transcripts : {len(block_dict)}\n"
          f"Number of nodes       : {len(nodes)}\n"
          f"--------------------------------")

    gfa_out = open(args.gfa, 'w')    
    print_log(f"Writing segments")
    for node in nodes:
        gfa_out.write(f"S\t{nodes[node]+1}\t{node}\tLN:i:{len(node)}\n")
    gfa_out.flush()

    print_log(f"Writing links")
    for transcript_name in block_dict:
        for seq in enumerate(block_dict[transcript_name].block_seqs):
            block_dict[transcript_name].path.append(nodes[seq[1]]+1)
        block_dict[transcript_name].get_links(gfa_out)
    gfa_out.flush()
    
    print_log(f"Writing paths")
    for transcript_name in block_dict:
        block_dict[transcript_name].write_path(gfa_out)
    gfa_out.flush()
    gfa_out.close()

    print_log(f"Creating clusters file")
    edges = pd.read_csv(f"{args.out}.edges.tsv", sep='\t')
    G = nx.from_pandas_edgelist(edges, source='query', target='target')
    components = nx.connected_components(G)
    output = open(f"{args.out}.clusters.tsv", 'w')
    # Output nodes of each component into separate files
    clustered_transcripts = []
    for i, component in enumerate(components):
        g_id = f'group{i}'
        for node in component:
            clustered_transcripts.append(node)
            output.write(f'{g_id}\t{node}\t{transcript_direction[node]}\t{len(seq_dict[node])}\tc\t{seq_dict[node].description}\n')
    
    for seq in seq_dict:
        if seq not in clustered_transcripts:
            i += 1
            g_id = f'group{i}'
            output.write(f'{g_id}\t{seq}\t+\t{len(seq_dict[seq])}\tu\t{seq_dict[seq].description}\n')
    output.close()

def index_dict(lst):
    """
    Index a list
    """
    return {v: k for k, v in enumerate(lst)}

def time_from_start():
    """
    Get the time from start
    """
    return f"{seconds_to_hhmmss(int(time.time() - start_time))}"

def export_path(args):
    """
    Export path from gfa file
    """
    parser = GFAParser(args.gfa)
    gfa = parser.parse()
    print(f"{gfa.summary()}")
    print_log(f"Extracting path: {args.path}")
    gfa.extract_path(args.path, args.out)

def sub_graph(args):
    """
    Extract sub graph from gfa file
    """
    parser = GFAParser(args.gfa)
    gfa = parser.parse()
    print(f"{gfa.summary()}")
    print_log(f"Extracting sub graph for node: {args.node}")
    sub_gfa = gfa.sub_graph(args.node)
    sub_gfa.write_to_file(args.out)

def extract_nodes(args):
    """
    Extract nodes from gfa file
    """
    parser = GFAParser(args.gfa)
    gfa = parser.parse()
    print(f"{gfa.summary()}")
    print_log(f"Extracting nodes: {args.nodes}")
    node_names = args.nodes.split(',')
    segments = [s for s in gfa.segments if s.name in node_names]
    with open(args.out, 'w') as file:
        for segment in segments:
            file.write(f">{segment.name}\n")
            file.write(f"{fold(segment.sequence, 60)}\n")

def aln_to_gfa(args):
    """
    Align fastq reads to gfa file
    """
    check_vg_installed()
    if args.fastq:
        fastq = args.fastq.split(',')
        if len(fastq) == 1:
            fastq.append(None)
    # check if index files exist
    if not os.path.isfile(f"{args.index}.xg") or not os.path.isfile(f"{args.index}.gcsa"):
        print(f"Index files not found for {args.index}. Please run 'wp2_test_gfa.py index -i {args.index} -w map'")
        sys.exit(1)

    if args.out.endswith('.gam'):
        # do aling with vg map
        print_log("Aligning reads to graph")
        if fastq[1]:
            run_cmd(f"vg map -x {args.index}.xg -g {args.index}.gcsa -t {args.threads} -f {fastq[0]} -f {fastq[1]} > {args.out}")
        else:
            run_cmd(f"vg map -x {args.index}.xg -g {args.index}.gcsa -t {args.threads} -f {fastq[0]} > {args.out}")   
        
        # # do align using vg giraffe, check if index files exist with .dist extension and .min extension.
        # print_log("Aligning reads to graph")
        # if fastq[1]:
        #     run_cmd(f"vg giraffe -x {args.index}.xg -t {args.threads} --fastq-in {fastq[0]} --fastq-in {fastq[1]} -o gaf > {args.out}")
        # else:
        #     run_cmd(f"vg giraffe -x {args.index}.xg -t {args.threads} --fastq-in {fastq[0]} -o gaf > {args.out}")
        if args.json:
            # convert gamp to json
            print_log("Converting gam to json")
            run_cmd(f"vg view --threads {args.threads} -aj {args.out} > {args.out}.json")

    elif args.out.endswith('.gamp'):
        # do align using vg mpmap
        print_log("Aligning reads to graph")
        if fastq[1]:
            run_cmd(f"vg mpmap --not-spliced --hit-max 100 -t {args.threads} -x {args.index}.xg -g {args.index}.gcsa -n DNA -f {fastq[0]} -f {fastq[1]} > {args.out}")
        else:
            run_cmd(f"vg mpmap --not-spliced --hit-max 100 -t {args.threads} -x {args.index}.xg -g {args.index}.gcsa -n DNA -f {fastq[0]} > {args.out}")
        
        if args.json:
            # convert gamp to json
            print_log("Converting gamp to json")
            run_cmd(f"vg view --threads {args.threads} -Kj {args.out} > {args.out}.json")
    else:
        print(f"Output file must be either .gam or .gamp")
        sys.exit(1)

def is_newer(file1, file2):
    """
    Check if file1 is newer than file2
    """
    ctime1 = os.path.getctime(file1)
    ctime2 = os.path.getctime(file2)
    return ctime1 > ctime2

def check_vg_installed():
    """
    check if tool vg is installed in PATH
    """
    try:
        subprocess.run(["vg"], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    except FileNotFoundError:
        print("Tool vg not found in PATH. Please install vg")
        sys.exit(1)

def run_cmd(cmd_line):
    """
    Run a command line
    """
    print_log(f"CMD: {cmd_line}")
    p = subprocess.Popen(cmd_line, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = p.communicate()
    if p.returncode != 0:
        raise Exception(err)
    return out.decode('utf-8'), err.decode('utf-8')

def fold(string, width):
    """
    Fold a string to a given width
    """
    return "\n".join([string[i:i+width] for i in range(0, len(string), width)])

def seconds_to_hhmmss(seconds):
    """
    Convert seconds to hh:mm:ss format
    """
    m, s = divmod(seconds, 60)
    h, m = divmod(m, 60)
    return f"{h:02d}:{m:02d}:{s:02d}"

def index_gfa(args):
    if not args.gfa.endswith('.gfa'):
        print(f"Input file must be a gfa file. '.gfa' extension is required")
        sys.exit(1)
    
    out = args.gfa.replace('.gfa', '')
    if not os.path.isfile(f"{out}.vg") or is_newer(f"{args.gfa}", f"{out}.vg"):
        print_log(f"Converting gfa to vg")
        run_cmd(f"vg convert --threads {args.threads} -g {args.gfa} -p > {out}.vg")
    else:
        print_log(f"Found vg file: {out}.vg")
        print_log(f"Skipping conversion")

    if not os.path.isfile(f"{out}.xg") or is_newer(f"{out}.vg", f"{out}.xg"):
        print_log(f"Indexing vg -> xg")
        run_cmd(f"vg index --threads {args.threads} -x {out}.xg {out}.vg")
    else:
        print_log(f"Found xg file: {out}.xg")
        print_log(f"Skipping indexing .xg")

    if not os.path.isfile(f"{out}.gcsa") or is_newer(f"{out}.xg", f"{out}.gcsa"):
        print_log(f"Indexing vg -> gcsa")
        run_cmd(f"vg prune --threads {args.threads} -k45 -M 32 -r {out}.vg > {out}.pruned.vg")
        run_cmd(f"vg index --threads {args.threads} -g {out}.gcsa -b ./tmp -Z 20 {out}.pruned.vg")
        run_cmd(f"rm {out}.pruned.vg")
    else:
        print_log(f"Found gcsa file: {out}.gcsa")
        print_log(f"Skipping indexing .gcsa")
    
    # if not os.path.isfile(f"{out}.dist") or is_newer(f"{out}.vg", f"{out}.dist"):
    #     print(f"Indexing vg -> dist")
    #     run_cmd(f"vg index --threads {args.threads} -j {out}.dist {out}.vg")

    # if not os.path.isfile(f"{out}.gbwt") or is_newer(f"{out}.xg", f"{out}.gbwt"):
    #     print(f"Indexing xg -> gbwt")
    #     run_cmd(f"vg gbwt --num-threads {args.threads} -x {out}.xg --index-paths -o {out}.gbwt")

    # if not os.path.isfile(f"{out}.min") or is_newer(f"{out}.vg", f"{out}.min"):
    #     print(f"Indexing vg -> min")
    #     run_cmd(f"vg minimizer --threads {args.threads} -g {out}.gbwt -d {out}.dist -o {out}.min {out}.vg")

def print_log(msg):
    """
    Print log messages
    """
    print(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] {msg}", file=sys.stderr)

def check_installed(tool):
    """
    Check if a tool is installed
    """
    try:
        subprocess.run([tool], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    except FileNotFoundError:
        print(f"Tool {tool} not found in PATH. Please install {tool}")
        sys.exit(1)

def interleave_gaf(gaf1, gaf2, gaf):
    """
    Interleave two gaf files
    """
    print_log(f"Interleaving gaf files")
    with open(gaf1, 'r') as f1, open(gaf2, 'r') as f2, open(gaf, 'w') as fout:
        for line1, line2 in zip(f1, f2):
            fout.write(line1)
            fout.write(line2)

def aln_to_gfa_ga(args):
    """
    Align fastq reads to gfa file using GraphAligner
    """
    check_installed('GraphAligner')
    print_log("Aligning reads to graph using GraphAligner")
    if len(args.fastq) == 1:
        args.fastq.append(None)
    elif len(args.fastq) > 2:
        print_log(f"ERROR: Maximum two fastq files are allowed")
        sys.exit(1)
    if args.fastq[1]:
        run_cmd(f"GraphAligner -g {args.gfa} -f {args.fastq[0]} -t {args.threads} -a {args.out}.1.gaf -x vg --multimap-score-fraction 1")
        run_cmd(f"GraphAligner -g {args.gfa} -f {args.fastq[1]} -t {args.threads} -a {args.out}.2.gaf -x vg  --multimap-score-fraction 1")
        # interleave_gaf(f"{args.out}.1.gaf", f"{args.out}.2.gaf", f"{args.out}.gaf")
    else:
        run_cmd(f"GraphAligner -g {args.gfa} -f {args.fastq[0]} -t {args.threads} -a {args.out}.gaf -x vg  --multimap-score-fraction 1")

def update_strand_info(transcript_direction, qName, tName, strand):
    """
    Update the strand information from the psl file
    """
    # Directionality of both transcripts already defined
    if((tName in transcript_direction) and (qName in transcript_direction)):
        return transcript_direction
    # Transcript directiionality already defined
    elif(tName in transcript_direction): 
        tdir = transcript_direction[tName]
        if(tdir == strand): #If both the transcript and strand same
            transcript_direction[qName] = '+'
        else:
            transcript_direction[qName] = '-'
    # Query transcript already defined
    elif(qName in transcript_direction): 
        qdir = transcript_direction[qName]
        if(qdir == strand): # If both the query and strand same
            transcript_direction[tName] = '+'
        else:
            transcript_direction[tName] = '-'
    # Neither of the sequences defined
    else:
        if(strand=='+'): # Both transcripts in the same direction
            transcript_direction[qName] = '+'
            transcript_direction[tName] = '+'

        else: # arbitrarily define transcript as postive and query as negative
            transcript_direction[tName] = '+'
            transcript_direction[qName] = '-'
    return transcript_direction

def get_sets(ilist, n):
    olist = []
    for i in range(len(ilist)-n+1):
        olist.append(tuple(ilist[i:i+n]))
    return olist


def counts_from_gaf(args):
    """
    Get counts from gaf file for each path
    """
    parser = GFAParser(args.gfa)
    gfa = parser.parse()
    print(f"{gfa.summary()}")
    segments = {s.name: s for s in gfa.segments}
    paths = {p.name: p for p in gfa.paths}
    tr_to_cltr = defaultdict(list)
    node_to_clusters = defaultdict(list)
    counts = defaultdict(lambda: defaultdict(int))

    print_log(f"Parsing clusters file: {args.clusters}")
    with open(args.clusters, 'r') as file:
        for line in file:
            fields = line.strip().split("\t")
            tr_to_cltr[fields[1]] = fields[0]
            if fields[0] not in counts:
                counts[fields[0]] = {'count': 0, 'ambigus': 0, 'lengths': []}
            counts[fields[0]]['lengths'].append(int(fields[3]))
            

    path_index = defaultdict(dict)
    for path in paths:
        if len(paths[path].segments) < 2:
            continue
        for i in range(2, 5):
            if len(paths[path].segments) < i:
                break
            for node_set in get_sets(paths[path].segments, i):
                path_index[node_set] = path
            
    
    for path in paths:    
        for node in paths[path].segments:
            if tr_to_cltr[path] not in node_to_clusters[node]:
                node_to_clusters[node].append(tr_to_cltr[path])
            
    for in_gaf in args.input:
        print_log(f"Extracting counts from gaf file: {in_gaf}")
        with open(in_gaf, 'r') as file:
            first_gaf = True
            for line in file:
                clusters = []
                gaf = GAF(line)
                if not first_gaf:
                    if last_gaf.read_name == gaf.read_name:
                        continue
                first_gaf = False
                nodes = gaf.nodes()
                if len(nodes[1]) == 1:
                    if nodes[1][0] in node_to_clusters:
                        clusters.extend(node_to_clusters[nodes[1][0]])
                    else:
                        print(f"Node {nodes[1][0]} not found in gfa file")
                elif tuple(nodes[1]) in path_index:
                    path = path_index[tuple(nodes[1])]
                    clusters.append(tr_to_cltr[path])
                if len(clusters) == 1:
                    counts[clusters[0]]['count'] += 1
                elif len(clusters) > 1:
                    for cluster in set(clusters):
                        counts[cluster]['ambigus'] += 1
                last_gaf = gaf

    print_log(f"Writing counts to file: {args.out}.counts")
    ofh = open(f"{args.out}.counts", 'w')
    ofh.write(f"cluster\ttr_count\tmean_length\tcount\tambigus\n")
    for cluster in counts:
        ofh.write(f"{cluster}\t{len(counts[cluster]['lengths'])}\t{int(sum(counts[cluster]['lengths'])/len(counts[cluster]['lengths']))}\t{counts[cluster]['count']}\t{counts[cluster]['ambigus']}\n")
            

def main():
    """
    Main function
    """    
    # test commands
    # wp2_test_gfa.py psl_to_gfa -i rna.psl -f rna.fna -o rna.2.gfa -c 0.8 -b 4
    # wp2_test_gfa.py index -i rna.2.gfa -t 20
    # wp2_test_gfa.py aln_to_gfa -i rna.2 -f rna.read1.fq,rna.read2.fq -o rna.2.gam -t 20 -j

    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest='command')

    # parse arguments for the "psl_to_gfa" command
    psl_to_gfa_parser = subparsers.add_parser('psl_to_gfa', help='Convert psl file to gfa file')
    psl_to_gfa_parser.add_argument('-i', '--psl', help='blast psl output file', required=True)
    psl_to_gfa_parser.add_argument('-f', '--fasta', help='fasta file', required=True)
    psl_to_gfa_parser.add_argument('-o', '--gfa', help='gfa file', required=True)
    psl_to_gfa_parser.add_argument('-p', '--outpsl', help='output filtered psl file')
    psl_to_gfa_parser.add_argument('-c', '--coverage', help='minimum coverage 0.0-1.0 [0.8]', default=0.8, type=float)
    psl_to_gfa_parser.add_argument('-b', '--blocks', help='maximum number of blocks [4]', default=4, type=int)

    # parse arguments for the "export_path" command
    export_path_parser = subparsers.add_parser('export_path', help='Export path from gfa file')
    export_path_parser.add_argument('-i', '--gfa', help='gfa file', required=True)
    export_path_parser.add_argument('-o', '--out', help='output fasta file', required=True)
    export_path_parser.add_argument('-p', '--path', help='path name', required=True)

    # parse arguments for the "sub_graph" command
    sub_graph_parser = subparsers.add_parser('sub_graph', help='Extract sub graph from gfa file')
    sub_graph_parser.add_argument('-i', '--gfa', help='gfa file', required=True)
    sub_graph_parser.add_argument('-n', '--node', help='node name', required=True)
    sub_graph_parser.add_argument('-o', '--out', help='output gfa file', required=True)

    # parse arguments for "extract_nodes" command
    extract_nodes_parser = subparsers.add_parser('extract_nodes', help='Extract nodes from gfa file')
    extract_nodes_parser.add_argument('-i', '--gfa', help='gfa file', required=True)
    extract_nodes_parser.add_argument('-o', '--out', help='output fasta file', required=True)
    extract_nodes_parser.add_argument('-n', '--nodes', help='node names (one or more comma seperated)', required=True)

    # parse arguments for "index" command
    index_parser = subparsers.add_parser('index', help='Index gfa file')
    index_parser.add_argument('-i', '--gfa', help='gfa file', required=True)
    index_parser.add_argument('-t', '--threads', help=f"number of threads [{nthreads}]", default=nthreads, type=int)
    
    # parse arguments for "aln_to_gfa" command
    aln_to_gfa_parser = subparsers.add_parser('aln_to_gfa', help='algin fastq reads to gfa file')
    aln_to_gfa_parser.add_argument('-i', '--index', help='index prefix of .xg and .gcsa', required=True)
    aln_to_gfa_parser.add_argument('-f', '--fastq', help='fastq file, if paired end comma seperated', required=True)
    aln_to_gfa_parser.add_argument('-o', '--out', help='output file. either .gam or .gamp', required=True)
    aln_to_gfa_parser.add_argument('-t', '--threads', help=f"number of threads [{nthreads}]", default=nthreads, type=int)
    aln_to_gfa_parser.add_argument('-j', '--json', help='output json file', action='store_true')

    # parse arguments for "aln_to_gfa_ga" command
    aln_to_gfa_ga_parser = subparsers.add_parser('aln_to_gfa_ga', help='algin fastq reads to gfa file using GraphAligner')
    aln_to_gfa_ga_parser.add_argument('-g', '--gfa', help='input GFA', required=True)
    aln_to_gfa_ga_parser.add_argument('-f', '--fastq', help='fastq file, if paired end comma seperated', nargs='+', required=True)
    aln_to_gfa_ga_parser.add_argument('-o', '--out', help='output file prefix. .gaf', required=True)
    aln_to_gfa_ga_parser.add_argument('-t', '--threads', help=f"number of threads [{nthreads}]", default=nthreads, type=int)

    # parse arguments for "counts_from_gaf" command
    counts_from_gaf_parser = subparsers.add_parser('counts_from_gaf', help='count reads per cluster from gaf file')
    counts_from_gaf_parser.add_argument('-i', '--input', help='GraphAligner output gaf file', nargs='+', required=True)
    counts_from_gaf_parser.add_argument('-g', '--gfa', help='graph gfa file', required=True)
    counts_from_gaf_parser.add_argument('-c', '--clusters', help='clusters file', required=True)
    counts_from_gaf_parser.add_argument('-o', '--out', help='output prefix', required=True)

    # parse arguments for "aln" command
    aln_parser = subparsers.add_parser('aln', help='align reads to graph')
    aln_parser.add_argument('-g', '--gfa', help='graph gfa file', required=True)
    aln_parser.add_argument('-f', '--fastq', help='fastq file, if paired end comma seperated', nargs='+', required=True)
    aln_parser.add_argument('-o', '--out', help='output file prefix. .gaf', required=True)
    aln_parser.add_argument('-t', '--threads', help=f"number of threads [{nthreads}]", default=nthreads, type=int)

    args = parser.parse_args(args=(sys.argv[1:] or ['--help']))

    if args.command == 'psl_to_gfa':
        psl_to_gfa(args)
    
    if args.command == 'export_path':
        export_path(args)

    if args.command == 'sub_graph':    
        sub_graph(args)
    
    if args.command == 'extract_nodes':
        extract_nodes(args)
    
    if args.command == 'index':
        index_gfa(args)

    if args.command == 'aln_to_gfa':
        aln_to_gfa(args)
    
    if args.command == 'aln_to_gfa_ga':
        aln_to_gfa_ga(args)

    if args.command == 'counts_from_gaf':
        counts_from_gaf(args)
    
    # if args.command is 'aln':
    #     aln(args)
    
    # end clock
    end_time = time.time()
    print_log(f"Time taken: {seconds_to_hhmmss(int(end_time - start_time))}")


if __name__ == "__main__":
    main()

# EOF