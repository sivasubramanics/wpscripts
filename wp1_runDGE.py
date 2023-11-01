#!/usr/bin/env python
import sys, os
from collections import defaultdict
import pandas as pd
import subprocess

# class object carrying sample information with its factors and replicates
class Sample():
    def __init__(self, name):
        self.name = name
        self.reps = []
        self.factors = defaultdict()
    
    def get_reps(self):
        return self.reps

    def set_reps(self, repname):
        self.reps.append(repname)
    
    def set_factor(self, factor_key, factor_value):
        if factor_key not in self.factors:
            self.factors[factor_key] = factor_value
        # else:
            # print(f"{self.name} already has a value for {factor_key}: {self.factors[factor_key]}")
    
    def get_factor_dict(self):
        return self.factors
        
    def get_factors_list(self):
        return list(self.factors.values())
    
    def get_factor(self, factor_key):
        return self.factors[factor_key]
    
    def __str__(self):
        return self.name

    def __repr__(self):
        return self.name
    
# process the conditions file and create a dictionary of sample objects
def get_conditions(in_conditions):
    with open(in_conditions, 'r') as fh:
        line_number = 0
        samples_dict = defaultdict()
        for line in fh:
            line = line.strip()
            data = line.split()
            line_number += 1
            if line_number == 1:
                head_data = data
                continue
            if data[1] not in samples_dict:
                samples_dict[data[1]] = Sample(data[1])
            samples_dict[data[1]].set_reps(data[0])
            for i in range(2, len(head_data)):
                samples_dict[data[1]].set_factor(head_data[i], data[i])
    return samples_dict

# def run_dge(sample_a, sample_b, count_matrix, dge_result, p_value, lfc, shfh):
#     r_script_path = os.path.abspath( os.path.dirname( __file__ ) ) + '/deseq2.R'
#     log_file = out_prefix + "_" + sample_a + "_vs_" + sample_b + "_DESeq2.log"
#     # fo = open(log_file, 'w')
#     # dge_run_log = subprocess.run(['Rscript', r_script_path, sample_a, sample_b, p_value, lfc], capture_output = True, text = True)
#     shfh.write(' '.join(['Rscript', r_script_path, sample_a, sample_b, count_matrix, dge_result, p_value, lfc, '>', log_file, '2>>', log_file]) + '\n')
#     # fo.write(f"{dge_run_log.stdout}\n")
#     # fo.write(f"{dge_run_log.stderr}\n")
#     # fo.close()
#     # if not dge_run_log.stderr:
#         # os.remove(dge_run_log)

# check if the dge is successful
def check_error(infile):
    with open(infile, 'r') as fh:
        for line in fh:
            line = line.strip()
            if 'error' in line or 'Error' in line or 'ERROR' in line:
                return True
    return False

# write shell command for running dge accoding to the method specified
def run_dge(sample_a, sample_b, count_matrix, dge_result, shfh, dge_mthd):
    r_script_path = os.path.abspath( os.path.dirname( __file__ ) ) + '/' + dge_mthd +'.R'
    log_file = out_prefix + "_" + sample_a + "_vs_" + sample_b + "_" + dge_mthd + ".log"
    shfh.write(' '.join(['Rscript', r_script_path, sample_a, sample_b, count_matrix, dge_result, '>', log_file, '2>>', log_file]) + '\n')
    
# get the contrast between two samples
def get_contrast(sample_a, sample_b):
    out_contrast = []
    keys_a = sample_a.keys()
    keys_b = sample_b.keys()
    if keys_a == keys_b:
        for key in keys_a:
            if sample_a[key] != sample_b[key]:
                out_contrast.append(key)
    return "; ".join(out_contrast)


# main function
if __name__ == "__main__":

    # command line arguments check
    if len(sys.argv) < 4:
        print(f"USAGE: {sys.argv[0]} <in_matrix> <conditions_file.txt> <output_prefix> <n_cpus> <deg_method>")
        exit(1)

    in_matrix = sys.argv[1] # count matrix
    in_conditions = sys.argv[2] # conditions file
    # p_value = sys.argv[3]
    # lfc = sys.argv[4]
    out_prefix = sys.argv[3] # output prefix
    n_cpus = int(sys.argv[4]) # number of cpus
    dge_method = sys.argv[5] # DESeq2 or edgeR
    out_file = out_prefix + '_' + dge_method +'_dge.tsv' # output file

    if dge_method == 'DESeq2':
        dge_mthd = 'deseq2'
    elif dge_method == 'edgeR':
        dge_mthd = 'edgeR'
    else:
        print(f"{dge_method} should be either DESeq2 or edgeR")
        exit(1)

    # read the conditions file. example as follows:
    # run     name    treatment       species time
    # IL12BB  IL12B   treated il      12hrs
    # IL12BC  IL12B   treated il      12hrs
    # IL12BD  IL12B   treated il      12hrs
    samples_dict = get_conditions(in_conditions)
    samples = list(samples_dict.keys())

    # read the count matrix
    count_df = pd.read_csv(in_matrix, sep = "\t")
    count_df.rename(columns = {count_df.columns[0]:'trans_id'}, inplace = True) 

    # choose the output file header based on the dge method
    ofh = open(out_file, 'w')
    if dge_method == "DESeq2":
        ofh.write("trans_id\tcontrast\tsampleA\tsampleB\tbaseMeanA\tbaseMeanB\tbaseMean\tlog2FoldChange\tlfcSE\tstat\tpvalue\tpadj\n")
    elif dge_method == "edgeR":
        ofh.write("trans_id\tcontrast\tsampleA\tsampleB\tlogFC\tlogCPM\tPValue\tFDR\n")
    ofh.close()

    # createa a shell script to run the dge commands
    shfh = open(out_prefix + '_dge.sh','w')
    for i in range(0, len(samples)):
        for j in range(i+1, len(samples)):
            sample_a = samples[i]
            sample_b = samples[j]
            if sample_a[:-1] != sample_b[:-1]:
                continue
            contrast = get_contrast(samples_dict[sample_a].get_factor_dict(), samples_dict[sample_b].get_factor_dict())
            print(f"{sample_a} vs {sample_b} : {contrast}")
            columns = ['trans_id']
            columns += samples_dict[sample_a].get_reps()
            columns += samples_dict[sample_b].get_reps()
            count_file = out_prefix + "_" + sample_a + "_vs_" + sample_b + "_count.matrix"
            dge_result = out_prefix + "_" + sample_a + "_vs_" + sample_b + "." + dge_mthd + ".DE_results.tsv"
            count_df.filter(columns).to_csv(count_file, index = False, sep = '\t')
            # run_dge(sample_a, sample_b, count_file, dge_result, p_value, lfc, shfh)
            run_dge(sample_a, sample_b, count_file, dge_result, shfh, dge_mthd)
    shfh.close()

    # read commands file and append them into a list
    with open(out_prefix + '_dge.sh') as file:
        commands = file.readlines()
        commands = [command.rstrip() for command in commands]

    # run the commands in parallel for the number of cpus provided
    for j in range(max(int(len(commands)/n_cpus)+1, 1)):
        procs = [subprocess.Popen(i, shell=True) for i in commands[j*n_cpus: min((j+1)*n_cpus, len(commands))] ]
        for p in procs:
            p.wait()

    # process the dge output and compile them according to the conditions and its factors
    for i in range(0, len(samples)):
        for j in range(i+1, len(samples)):
            sample_a = samples[i]
            sample_b = samples[j]
            if sample_a[:-1] != sample_b[:-1]:
                continue
            contrast = get_contrast(samples_dict[sample_a].get_factor_dict(), samples_dict[sample_b].get_factor_dict())
            columns = ['trans_id']
            columns += samples_dict[sample_a].get_reps()
            columns += samples_dict[sample_b].get_reps()
            log_file = out_prefix + "_" + sample_a + "_vs_" + sample_b + "_" + dge_mthd + ".log"
            if check_error(log_file):
                print(f"Error in running {dge_method}.R")
                exit(1)
            else:
                os.remove(log_file)
            count_file = out_prefix + "_" + sample_a + "_vs_" + sample_b + "_count.matrix"
            dge_result = out_prefix + "_" + sample_a + "_vs_" + sample_b + "." + dge_mthd + ".DE_results.tsv"
            result_df = pd.read_csv(dge_result, sep = "\t")
            result_df.rename(columns = {result_df.columns[0]:'trans_id'}, inplace = True)
            result_df.insert(1, 'contrast', contrast)
            result_df.to_csv(out_file, mode = 'a', header = False, index = False, sep = '\t')
            # os.remove(dge_result)
            os.remove(count_file)
            
            
