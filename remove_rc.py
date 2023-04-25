import os
import csv
import sys
import argparse
from pathlib import Path
import re
from Bio.Seq import Seq
from Bio import SeqIO

def check_args(args):
    if args.file == args.output:
        raise Exception("Input and output file cannot be the same.")
    return

def check_file_exists(output_filepath):
    file = Path(output_filepath)
    if file.is_file():
        reply = input(f"File {output_filepath} already exists. Overwrite? [y/n] ")
        if reply == "y":
            return
        else:
            print("Command terminated.")
            exit()

# def get_parameters(input_filepath):
#     params = {}
#     with open(input_filepath, "r") as f:
#         for line in f:
#             if line.startswith(";"):
#                 match = re.search("(L|k|bases|prevented_patterns|rcfree|nucleotide_composition)=(\d+|\w+|\"\w+\"|.+)", line)
#                 if match is not None:
#                     param_key = match.group(1)
#                     param_value = match.group(2)
#                     params[param_key] = param_value
#             elif line.startswith(">"):
#                 break
#             else:
#                 continue
#     return params

def add_subseq(seq_id, subseq, subseqs_dict, is_implied):
    if subseqs_dict.get(subseq) is None:
        subseqs_dict[subseq] = seq_id
    else:
        if is_implied:
            return subseqs_dict
        else:
            raise Exception("Subsequence already exists in another sequence. Sequences in input FASTA file not orthogonal, please regenerate sequences with generate_seqs.py")
    return subseqs_dict

def get_rc_graph(params):
    k = int(params.get("k"))
    subseqs_dict = {}
    rc_graph_dict = {}
    count = 0
    input_fasta = SeqIO.parse(open(params.get("input_filepath")),'fasta')

    for line in input_fasta:
        count += 1
        seq_id, seq_sequence = line.id, str(line.seq)
        rc_graph_dict[seq_id] = set()
        L = len(seq_sequence)
        for i in range(L-k+1):
            subseq = seq_sequence[i: i+k]
            subseqs_dict = add_subseq(seq_id, subseq, subseqs_dict, params.get("implied"))
            subseq_rc = str(Seq(subseq).reverse_complement())
            subseq_rc_seq_id = subseqs_dict.get(subseq_rc)
            # if RC is already a subseq encountered before (i.e. 2 sequences are linked by RC subsequences)
            if subseq_rc_seq_id is not None:
                rc_graph_dict[seq_id].add(subseq_rc_seq_id)
                rc_graph_dict[subseq_rc_seq_id].add(seq_id)
    params["size_before"] = count
    return rc_graph_dict

def trash_neighbours(seq_id, rc_graph, grouped_seqs):
        for neighbour in rc_graph.get(seq_id):
            grouped_seqs[neighbour] = "trashbin"

# def write_comments(params, output_filepath, exclude=list()):
#     for param_key, param_value in params.items():
#         if param_key in exclude:
#             continue
#         else:
#             output_filepath.write(f";{param_key}={param_value}\n")
#     output_filepath.write("\n")

def generate_rcfree_file(params, rcfree_seqs):
    input_fasta = SeqIO.parse(open(params.get("input_filepath")),'fasta')
    count = 0
    with open(params.get("output_filepath"), "w") as output_filepath:
        # write_comments(params, output_filepath, exclude=["input_filepath", "output_filepath", "size_before"])
        for line in input_fasta:
            seq_id, seq_sequence = line.id, str(line.seq)
            if seq_id in rcfree_seqs.keys():
                count += 1
                output_filepath.write(f">{seq_id}\n{seq_sequence}\n") 
    print(f"""Reversed complements have been removed! \nSize of library before RC removal: {params.get('size_before')} \nSize of library after RC removal: {count}""")       

def remove_rc(params):
    grouped_seqs = {}
    rc_graph = get_rc_graph(params)

    # run through seq list from least connected to most connected
    for seq_id in sorted(rc_graph, key=lambda k: len(rc_graph[k]), reverse=False):
        seq_id_group = grouped_seqs.get(seq_id)
        if seq_id_group is None:
            grouped_seqs[seq_id] = "keeper"
            trash_neighbours(seq_id, rc_graph, grouped_seqs)
        elif seq_id_group == "keeper":
            trash_neighbours(seq_id, rc_graph, grouped_seqs)
        elif seq_id_group == "trashbin":
            continue
        else:
            raise Exception(f"Sequence status of sequence {seq_id} is {seq_id_group}. It can only be 'keeper', 'trashbin', or None.")
        
    rcfree_seqs = {seq:True for seq, group in grouped_seqs.items() if group == "keeper"}
    generate_rcfree_file(params, rcfree_seqs)

def main():
    user_input = argparse.ArgumentParser()
    required_arguments = user_input.add_argument_group('required arguments')
    optional_arguments = user_input.add_argument_group('optional arguments')
    required_arguments.add_argument('-f', '--file', action='store', default=None, type=str, required=True,
                                    help='Path to input FASTA file')
    required_arguments.add_argument('-o', '--output', action='store', required=True,
                                    default=None, type=str, 
                                    help="""Specify the output FASTA filepath (e.g. my_sequeces_filtered.fa)
                                    """)
    required_arguments.add_argument('-k', '--k', action='store', required=True,
                                    default=None, type=int, 
                                    help="""Specify the value of k, where all reverse complements of any subsequence of length k will be removed (e.g. 13)
                                    """)
    optional_arguments.add_argument('-i', '--implied', action='store_true',
                                    help="""Specify as True if using for implied barcodes - if True, it skips the orthogonality check.""")
    args = user_input.parse_args()
    check_args(args)
    check_file_exists(args.output)
    # params = get_parameters(args.file)
    params = {}
    params["rcfree"] = True
    params["input_filepath"] = args.file
    params["output_filepath"] = args.output
    params["k"] = args.k
    params["implied"] = args.implied
    remove_rc(params)

if __name__ == "__main__":
    main()
    