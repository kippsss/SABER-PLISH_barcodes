import argparse
from pathlib import Path
import math
import random
import csv
from Bio.Seq import Seq
from Bio import SeqIO

def check_args(args):
    if not args.file.endswith(".fa") and not args.file.endswith(".fasta"):
        raise Exception("Input file must be a FASTA file (.fa).")
    if not args.output.endswith(".tsv"):
        raise Exception("Output file must be a TSV file (.tsv).")

def check_file_exists(output_filepath):
    file = Path(output_filepath)
    if file.is_file():
        reply = input(f"File {output_filepath} already exists. Overwrite? [y/n] ")
        if reply == "y":
            return
        else:
            print("Command terminated.")
            exit()

def get_complexity(params):
    input_fasta = SeqIO.parse(open(params.get("input_filepath")),'fasta')
    with open(params.get("output_filepath")+"copy", "w") as output_filepath:
        output_writer = csv.writer(output_filepath, delimiter="\t")
        output_writer.writerow(["BARCODE_ID", "BARCODE_SEQ", "COMPLEXITY"])
        for line in input_fasta:
            id, seq = line.id, str(line.seq)
            length = len(seq)
            complexity_list = []
            for word_size in range(1, length):
                denominator = min(length-word_size+1, 4 ** word_size)
                words = {}
                for i in range(length-word_size+1):
                    word = seq[i: i+word_size]
                    if word not in words.keys():
                        words[word] = 1
                numerator = len(words)
                complexity_list.append(numerator/denominator)
            complexity = math.prod(complexity_list)
            output_writer.writerow([id, seq, complexity])

def main():
    user_input = argparse.ArgumentParser()
    required_arguments = user_input.add_argument_group('required arguments')
    required_arguments.add_argument('-f', '--file', action='store', default=None, type=str, required=True,
                                    help='Path to input FASTA file')
    required_arguments.add_argument('-o', '--output', action='store', required=True,
                                    default=None, type=str, 
                                    help="""Specify the output TSV filepath (e.g. 42mer_k14_complexity.tsv)
                                    """)
    args = user_input.parse_args()
    check_args(args)
    params = {}
    params["input_filepath"] = args.file
    params["output_filepath"] = args.output
    get_complexity(params)

if __name__ == '__main__':
    main()