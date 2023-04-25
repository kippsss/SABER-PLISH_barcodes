import argparse
from pathlib import Path
import math
import random
import csv
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqUtils import MeltingTemp as mt

def check_args(args):
    if not args.file.endswith(".fa") and not args.file.endswith(".fasta"):
        raise Exception("Input file must be a FASTA file (.fa).")
    if not args.output.endswith(".tsv"):
        raise Exception("Output file must be a TSV file (.tsv).")
    
def get_Tm(params):
    input_fasta = SeqIO.parse(open(params.get("input_filepath")),'fasta')
    with open(params.get("output_filepath"), "w") as output_filepath:
        output_writer = csv.writer(output_filepath, delimiter="\t")
        output_writer.writerow(["BARCODE_ID", "BARCODE_SEQ", "TM_ROUND", "TM"])
        count = 0
        for line in input_fasta:
            id, seq = line.id, str(line.seq)
            Tm = mt.Tm_NN(Seq(seq), nn_table=mt.DNA_NN2, Na=390)
            Tm_round = round(Tm ,2)
            output_writer.writerow([id, seq, Tm_round, Tm])
            count += 1
    print(f"Calculation of melting temperature for {count} sequences completed. ")


def main():
    user_input = argparse.ArgumentParser()
    required_arguments = user_input.add_argument_group('required arguments')
    required_arguments.add_argument('-f', '--file', action='store', default=None, type=str, required=True,
                                    help='Path to input FASTA file')
    required_arguments.add_argument('-o', '--output', action='store', required=True,
                                    default=None, type=str, 
                                    help="""Specify the output TSV filepath (e.g. 42mer_k14_Tm.tsv)
                                    """)
    args = user_input.parse_args()
    check_args(args)
    params = {}
    params["input_filepath"] = args.file
    params["output_filepath"] = args.output
    get_Tm(params)

if __name__ == '__main__':
    main()