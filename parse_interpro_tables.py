#!/usr/bin/env python3
# auther: siva.selvanayagam@wur.nl
# description: This script parses the interproscan tsv file and extracts the domain information for a given database
# usage: python3 parse_interpro_tables.py <interpro_tsv_file> <database>
# example: python3 parse_interpro_tables.py interproscan.tsv Pfam
# date: 03-05-2023

import os
import pathlib
import shutil
import sys
import urllib.request
from collections import defaultdict


class GO:
    """
    This class is used to store the GO information
    """

    def __init__(self, go_id):
        self.go_id = go_id
        self.go_term = ""
        self.go_namespace = ""
        self.go_definition = ""

    def set_go_term(self, go_term):
        self.go_term = go_term

    def set_go_namespace(self, go_namespace):
        self.go_namespace = go_namespace

    def set_go_definition(self, go_definition):
        self.go_definition = go_definition

    def get_go_term(self):
        return self.go_term

    def get_go_namespace(self):
        return self.go_namespace

    def get_go_definition(self):
        return self.go_definition

    def to_string(self):
        return self.go_id + "\t" + self.go_term + "\t" + self.go_definition

    def __str__(self):
        return self.go_id

    def __repr__(self):
        return self.go_id


def parse_obo_file(obo_file):
    """This function parses the obo file and returns a dictionary with GO id as key and GO term as value
        :param obo_file: obo file
        :return: dictionary with GO id as key and GO term as value
    """
    print(f"Parsing obo file {obo_file}")
    go_dict = defaultdict()
    with open(obo_file, 'r') as fh:
        for line in fh:
            if line.startswith("id:"):
                go_id = line.rstrip().split(": ")[1]
                go_dict[go_id] = GO(go_id)
            if line.startswith("name:"):
                go_term = line.rstrip().split(": ")[1]
                go_dict[go_id].set_go_term(go_term)
            if line.startswith("namespace:"):
                go_namespace = line.rstrip().split(": ")[1]
                go_dict[go_id].set_go_namespace(go_namespace)
            if line.startswith("def:"):
                go_definition = line.rstrip().split(": ")[1]
                go_dict[go_id].set_go_definition(go_definition)
    return go_dict


def download_file(url, file_name):
    """This function downloads the file from the url and saves it to the file_name
        :param url: url of the file to be downloaded
        :param file_name: name of the file to be saved
        :return: None
    """
    print(f"Downloading file {url}")
    with urllib.request.urlopen(url) as response, open(file_name, 'wb') as out_file:
        shutil.copyfileobj(response, out_file)


DBs = ["Coils", "Gene3D", "MobiDBLite", "PANTHER", "Pfam", "SUPERFAMILY", "TIGRFAM"]

if len(sys.argv) != 3:
    print("Usage: python3 parse_interpro_tables.py <interpro_tsv_file> <database>")
    sys.exit(1)

interpro_tsv_file = sys.argv[1]
database = sys.argv[2]
databases = []
if "," in database:
    databases = database.split(',')
    for database in databases:
        if database not in DBs:
            print("Database not recognized. Please choose from the following: " + str(DBs))
            sys.exit(1)
else:
    if database not in DBs:
        print("Database not recognized. Please choose from the following: " + str(DBs))
        sys.exit(1)
    databases.append(database)

go_db_file = pathlib.Path(sys.argv[0]).parent.absolute() / "data" / "go-basic.obo"
if not os.path.exists(go_db_file):
    os.makedirs(os.path.dirname(go_db_file), exist_ok=True)
    download_file("http://current.geneontology.org/ontology/go-basic.obo", go_db_file)

go_db_dict = parse_obo_file(go_db_file)
go_dict = defaultdict()
ann_dict = defaultdict()
with open(interpro_tsv_file, 'r') as fh:
    print(f"Parsing interproscan tsv file {interpro_tsv_file}")
    for line in fh:
        data = line.rstrip().split("\t")
        if data[3] not in databases:
            continue

        if data[3] not in ann_dict:
            ann_dict[data[3]] = defaultdict()
        if data[0] not in ann_dict[data[3]]:
            ann_dict[data[3]][data[0]] = []
        if data[5] not in ann_dict[data[3]][data[0]]:
            ann_dict[data[3]][data[0]].append(data[5])
        if data[0] not in go_dict:
            go_dict[data[0]] = []
        if len(data) > 13:
            gos = data[13].split("|")
            for go in gos:
                if go == "-":
                    continue
                if go not in go_dict[data[0]]:
                    go_dict[data[0]].append(go)

print(f"output files:\n")
for database in databases:
    out_prefix = interpro_tsv_file.replace(".fasta.tsv", "", 1) + "." + database + "."
    print(f"{database} ann\t- {out_prefix}annotation.tsv")
    fwd = open(out_prefix + "annotation.tsv", 'w')
    for key, ann in ann_dict[database].items():
        fwd.write(key + "\t" + "; ".join(ann) + "\n")
    fwd.close()


out_prefix = interpro_tsv_file.replace(".fasta.tsv", "", 1) + "."
print(f"GO id file\t- {out_prefix}go_annotation.tsv\n"
        f"GO MF ann file\t- {out_prefix}go_mf.tsv\n"
        f"GO BP ann file\t- {out_prefix}go_bp.tsv\n"
        f"GO CC ann file\t- {out_prefix}go_cc.tsv\n")

fwa = open(out_prefix + "go_annotation.tsv", 'w')
fw_mf = open(out_prefix + "go_mf.tsv", 'w')
fw_mf.write("# id\tgo_id\tgo_term\tgo_definition\n")
fw_bp = open(out_prefix + "go_bp.tsv", 'w')
fw_bp.write("# id\tgo_id\tgo_term\tgo_definition\n")
fw_cc = open(out_prefix + "go_cc.tsv", 'w')
fw_cc.write("# id\tgo_id\tgo_term\tgo_definition\n")
for key, gos in go_dict.items():
    if len(gos) > 0:
        fwa.write(key + "\t" + ";".join(gos) + "\n")
    for go in gos:
        if go_db_dict[go].get_go_namespace() == "molecular_function":
            fw_mf.write(key + "\t" + go_db_dict[go].to_string() + "\n")
        elif go_db_dict[go].get_go_namespace() == "biological_process":
            fw_bp.write(key + "\t" + go_db_dict[go].to_string() + "\n")
        elif go_db_dict[go].get_go_namespace() == "cellular_component":
            fw_cc.write(key + "\t" + go_db_dict[go].to_string() + "\n")
        else:
            print("GO namespace not recognized: " + go_db_dict[go].get_go_namespace())
            sys.exit(1)
fwa.close()
fw_mf.close()
fw_bp.close()
fw_cc.close()

fw_background = open(out_prefix + "go_background.tsv", 'w')
fw_background.write("\n".join(go_dict.keys()) + "\n")
