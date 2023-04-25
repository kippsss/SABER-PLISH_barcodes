import argparse
from pathlib import Path
import math
import random
import csv
from Bio.Seq import Seq
from Bio import SeqIO

def check_args(args):
    if args.file == args.output:
        raise Exception("Input and output file cannot be the same.")
    if not args.file.endswith(".fa") and not args.file.endswith(".fasta"):
        raise Exception("Input file must be a FASTA file (.fa).")
    
def check_file_exists(output_filepath):
    file = Path(output_filepath)
    if file.is_file():
        reply = input(f"File {output_filepath} already exists. Overwrite? [y/n] ")
        if reply == "y":
            return
        else:
            print("Command terminated.")
            exit()

def get_bc_combi(params):
    input_fasta = SeqIO.parse(open(params.get("input_filepath")),'fasta')
    
    with open(params.get("tsv_input_filepath"), "r") as tsv_filepath:
        with open(params.get("output_filepath"), "w") as output_filepath:
            output_writer = csv.writer(output_filepath, delimiter="\t")
            output_writer.writerow(["NB_BARCODE_ID", "NB_BARCODE_SEQ", "PROBE_BARCODE_ID", "PROBE_BARCODE_SEQ", "IMPLIED_BARCODE_ID", "IMPLIED_BARCODE_SEQ"])
            
            count = 0
            implied_bc_dict = {}
            for line in input_fasta:
                id, seq = line.id, str(line.seq)
                if implied_bc_dict.get(id) is None:
                    implied_bc_dict[id] = seq
                # if implied_bc_dict.get(seq) is None:
                #     implied_bc_dict[seq] = 1
                else:
                    raise Exception(f"Implied barcode duplicate encountered.\nImplied barcode ID:{id}\nImplied barcode ID:{seq}")
        
            all_bc_pre = csv.reader(tsv_filepath, delimiter="\t")
            past_header = False
            for row in all_bc_pre:
                if not past_header:
                    past_header = True
                    continue
                implied_bc_id = row[4]
                implied_bc_seq = row[5]
                if implied_bc_dict.get(implied_bc_id) == implied_bc_seq:
                    count += 1
                    output_writer.writerow(row)

    print(f"Final number of nanobody-probe-implied barcode combinations: {count}")

def main():
    user_input = argparse.ArgumentParser()
    required_arguments = user_input.add_argument_group('required arguments')
    required_arguments.add_argument('-f', '--file', action='store', default=None, type=str, required=True,
                                    help='Path to input FASTA file')
    required_arguments.add_argument('-t', '--tsv_input', action='store', required=True,
                                    default=None, type=str, 
                                    help="""Specify the input TSV filepath (e.g. all_bc.tsv)
                                    """)
    required_arguments.add_argument('-o', '--output', action='store', required=True,
                                    default=None, type=str, 
                                    help="""Specify the output FASTA filepath (e.g. my_sequeces_filtered.fa)
                                    """)
    args = user_input.parse_args()
    check_args(args)
    check_file_exists(args.output)
    params = {}
    params["input_filepath"] = args.file
    params["tsv_input_filepath"] = args.tsv_input
    params["output_filepath"] = args.output
    get_bc_combi(params)

if __name__ == '__main__':
    main()