import argparse
from pathlib import Path
import math
import random
import csv
from Bio.Seq import Seq
import os

def check_dir_exists(output_filepath):
    dir_name = os.path.dirname(output_filepath)
    if not os.path.exists(dir_name):
        os.makedirs(dir_name)
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
    else:
        check_dir_exists(output_filepath)
    return

def check_even(nb_seq, probe_seq):
    nb_length = len(nb_seq)
    probe_length = len(probe_seq)
    if nb_length != probe_length:
        print(f"Length of nanobody barcode and probe barcode not the same.\nNanobody barcode length: {nb_length}\nProbe barcode length: {probe_length}")
        reply = input(f"Continue anyway? [y/n] ")
        if reply == "n":
            print("Command terminated.")
            exit()
    if nb_length % 2 != 0 or probe_length % 2 != 0:
        raise Exception(f"Both nanobody barcode and probe barcode lengths must be an even number.\nNanobody barcode length: {nb_length}\nProbe barcode length: {probe_length}")
    return nb_length, probe_length

def get_implied_simple(params):
    nb_seq = params.get("nb_seq")
    probe_seq = params.get("probe_seq")
    nb_length, probe_length = check_even(nb_seq, probe_seq)
    nb_half_length = int(nb_length / 2)
    nb_half = nb_seq[nb_half_length:]
    probe_half_length = int(probe_length / 2)
    probe_half = probe_seq[:probe_half_length]
    
    implied_seq_rc = nb_half + probe_half
    implied_seq = str(Seq(implied_seq_rc).reverse_complement())

    with open(params.get("combi_output"), "w") as combi_output_filepath:
        combi_output_filepath.write(f">nb-seq\n{nb_seq}\n")
        combi_output_filepath.write(f">probe-seq\n{probe_seq}\n")
        combi_output_filepath.write(f">implied-seq\n{implied_seq}\n")

def main():
    user_input = argparse.ArgumentParser()
    required_arguments = user_input.add_argument_group('required arguments')
    required_arguments.add_argument('-n', '--nb_seq', action='store', required=True,
                                    default=None, type=str, 
                                    help="""Specify the nanobody barcode sequence (e.g. ATCGGCAATT)
                                    """)
    required_arguments.add_argument('-p', '--probe_seq', action='store', required=True,
                                    default=None, type=str, 
                                    help="""Specify the probe barcode sequence (e.g. ACTTGCTCCG)
                                    """)
    required_arguments.add_argument('-o', '--combi_output', action='store', required=True,
                                    default=None, type=str, 
                                    help="""Specify the output FASTA file that will contain (1) nanobody barcode sequence, (2) probe barcode sequence, (3) implied barcode sequence
                                    """)
    
    args = user_input.parse_args()
    params = {}
    params["nb_seq"] = args.nb_seq
    params["probe_seq"] = args.probe_seq
    params["combi_output"] = args.combi_output
    check_file_exists(args.combi_output)
    get_implied_simple(params)

if __name__ == '__main__':
    main()