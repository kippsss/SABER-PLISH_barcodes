import os
import csv
import sys
import re
import argparse
from Bio.Seq import Seq
from Bio import SeqIO
from pathlib import Path

def check_args(args):
    if args.file == args.output:
        raise Exception("Input and output file cannot be the same.")
    arg_key_val = {"A": args.A, "T": args.T, "C": args.C, "G": args.G}
    if args.A is None and args.T is None and args.C is None and args.G is None:
        raise Exception(f"""No base content boundaries have been set.
        Please provide an upper and lower boundary for at least 1 base.
        Base content boundaries: {arg_key_val}""")
    for base, content_bounds in arg_key_val.items():
        if content_bounds is not None:
            if len(content_bounds) != 2:
                raise Exception(f"""Lower and upper boundaries need to be set for {base}.
                Exactly 2 inputs required.
                Current lower and upper boundaries for {base}: {content_bounds}.
                Required: [<lower>, <upper>]""")
            if content_bounds[0] > content_bounds[1]:
                raise Exception(f"""Upper boundary needs to be >= lower boundary for {base}.
                Current lower and upper boundaries for {base}: {content_bounds}.
                Required: <upper> >= <lower>""")
        else:
            continue

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
#                 match = re.search("(L|k|bases|prevented_patterns|rcfree)=(\d+|\w+|\"\w+\"|.+)", line)
#                 if match is not None:
#                     param_key = match.group(1)
#                     param_value = match.group(2)
#                     params[param_key] = param_value
#             elif line.startswith(">"):
#                 break
#             else:
#                 continue
#     return params

# def write_comments(params, output_filepath, exclude=list()):
#     for param_key, param_value in params.items():
#         if param_key in exclude:
#             continue
#         else:
#             output_filepath.write(f";{param_key}={param_value}\n")
#     output_filepath.write("\n")

def within_bounds(pct, bounds):
    if bounds is None or (pct >= bounds[0] and pct <= bounds[1]):
        return True
    else:
        return False
    
def filter_composition(params):
    input_fasta = SeqIO.parse(open(params.get("input_filepath")),'fasta')

    with open(params.get("output_filepath"), "w") as output_filepath:
        # write_comments(params, output_filepath, exclude=["input_filepath", "output_filepath"])

        input_count = 0
        count = 0
        for line in input_fasta:
            input_count += 1
            id, seq = line.id, str(line.seq)
            length = len(seq)
            is_within = True
            for nucleotide, bounds in params.get("nucleotide_composition").items():
                nucleotide_pct = (seq.count(nucleotide) / length) * 100
                is_within = within_bounds(nucleotide_pct, bounds)
                if not is_within:
                    break
            if is_within:
                count += 1
                output_filepath.write(f">{id}\n{seq}\n")
        print(f"Base content filtering completed! \nSize of library before filtering nucleotide composition: {input_count} \nSize of library after filtering nucleotide composition: {count}")
        
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
    optional_arguments.add_argument('-a', '--A', action='store', nargs='+',
                                    default=None, type=float, 
                                    help="""Range of A content to filter. Lower and upper bounds inclusive, separated by '-'
                                    e.g. 24.0-26.0
                                    """)
    optional_arguments.add_argument('-t', '--T', action='store', nargs='+',
                                    default=None, type=float, 
                                    help="""Range of T content to filter. Lower and upper bounds inclusive, separated by '-'
                                    e.g. 24.0-26.0
                                    """)
    optional_arguments.add_argument('-c', '--C', action='store', nargs='+',
                                    default=None, type=float, 
                                    help="""Range of C content to filter. Lower and upper bounds inclusive, separated by '-'
                                    e.g. 24.0-26.0
                                    """)
    optional_arguments.add_argument('-g', '--G', action='store', nargs='+',
                                    default=None, type=float, 
                                    help="""Range of G content to filter. Lower and upper bounds inclusive, separated by '-'
                                    e.g. 24.0-26.0
                                    """)
    
    args = user_input.parse_args()
    check_args(args)
    check_file_exists(args.output)
    # params = get_parameters(args.file)
    params = {}
    params["input_filepath"] = args.file
    params["output_filepath"] = args.output
    params["nucleotide_composition"] = {
        "A": args.A,
        "T": args.T,
        "C": args.C,
        "G": args.G,
    }
    filter_composition(params)

if __name__ == '__main__':
    main()