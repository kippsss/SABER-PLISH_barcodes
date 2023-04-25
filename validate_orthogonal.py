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
    
def validate_orthogonal(params):
    input_fasta = SeqIO.parse(open(params.get("input_filepath")),'fasta')
    
    k = params.get("k")
    subsequence_dict = {}
    count = 0
    for line in input_fasta:
        id, seq = line.id, str(line.seq)
        L = len(seq)
        for i in range(L-k+1):
            subsequence = seq[i: i+k]
            if subsequence_dict.get(subsequence) is None:
                subsequence_dict[subsequence] = id
                count += 1
            else:
                print("Subsequence present in more than 1 sequence.")
                print(f"Subsequence: {subsequence}")
                print(f"First sequence ID: {subsequence_dict.get(subsequence)}")
                print(f"Second sequence ID: {id}")
                raise Exception()
    print(f"Library of size {count} is orthogonal to the degree k = {k}")

def main():
    user_input = argparse.ArgumentParser()
    required_arguments = user_input.add_argument_group('required arguments')
    required_arguments.add_argument('-f', '--file', action='store', default=None, type=str, required=True,
                                    help='Path to input FASTA file')
    required_arguments.add_argument('-k', '--k', action='store', default=None, type=int, required=True,
                                    help='Specify k, the degree of orthogonality to check.')
    args = user_input.parse_args()
    check_args(args)
    params = {}
    params["input_filepath"] = args.file
    params["k"] = args.k
    validate_orthogonal(params)

if __name__ == '__main__':
    main()