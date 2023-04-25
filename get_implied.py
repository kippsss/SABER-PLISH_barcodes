import argparse
from pathlib import Path
import math
import random
from Bio.Seq import Seq
from Bio import SeqIO

def check_args(args):
    if args.nb_input == args.implied_output or args.probe_input == args.implied_output:
        raise Exception("Input and output files cannot be the same.")
    if not args.nb_input.endswith(".fa") and not args.nb_input.endswith(".fasta"):
        raise Exception("Nanobody barcode input file must be a FASTA file (.fa).")
    if not args.probe_input.endswith(".fa") and not args.probe_input.endswith(".fasta"):
        raise Exception("Probe barcode input file must be a FASTA file (.fa).")
    if not args.implied_output.endswith(".fa") and not args.implied_output.endswith(".fasta"):
        raise Exception("Implied barcode output file must be a FASTA file (.fa).")
    # if not args.implied_output.endswith(".tsv"):
    #     raise Exception("Nanobody-Probe-Implied barcode combi output file must be a TSV file (.tsv).")
    
def check_file_exists(output_filepath):
    file = Path(output_filepath)
    if file.is_file():
        reply = input(f"File {output_filepath} already exists. Overwrite? [y/n] ")
        if reply == "y":
            return
        else:
            print("Command terminated.")
            exit()

def check_length(nb_bc_seq, probe_bc_seq):
    nb_bc_length = len(nb_bc_seq)
    probe_bc_length = len(probe_bc_seq)
    if nb_bc_length != probe_bc_length:
        raise Exception(f"Length of nanobody barcode and probe barcode not the same.\nNanobody barcode length: {nb_bc_length}\nProbe barcode length: {probe_bc_length}")
    if nb_bc_length % 2 != 0 or probe_bc_length % 2 != 0:
        raise Exception(f"Both nanobody barcode and probe barcode lengths must be an even number.\nNanobody barcode length: {nb_bc_length}\nProbe barcode length: {probe_bc_length}")
    return nb_bc_length

def get_implied(params):

    with open(params.get("implied_output_filepath"), "w") as implied_output_filepath:
        implied_output_count = 0

        for nb_line in SeqIO.parse(open(params.get("nb_input_filepath")),'fasta'):
            for probe_line in SeqIO.parse(open(params.get("probe_input_filepath")),'fasta'):
                nb_id, nb_seq = nb_line.id, str(nb_line.seq)
                probe_id, probe_seq = probe_line.id, str(probe_line.seq)
                implied_id = f"{nb_id}_{probe_id}"
                length = check_length(nb_seq, probe_seq)
                half_length = int(length / 2)
                nb_half = nb_seq[half_length:]
                probe_half = probe_seq[:half_length]
                
                implied_seq_rc = nb_half + probe_half
                implied_seq = str(Seq(implied_seq_rc).reverse_complement())

                implied_output_filepath.write(f">{implied_id}\n{implied_seq}\n")
                implied_output_count += 1

    print(f"Number of Implied barcodes derived: {implied_output_count}")

def main():
    user_input = argparse.ArgumentParser()
    required_arguments = user_input.add_argument_group('required arguments')
    required_arguments.add_argument('-n', '--nb_input', action='store', required=True,
                                    default=None, type=str, 
                                    help="""Specify the nanobody barcode input FASTA filepath (e.g. 42mer_k14_nb.fa)
                                    """)
    required_arguments.add_argument('-p', '--probe_input', action='store', required=True,
                                    default=None, type=str, 
                                    help="""Specify the ISH probe barcode input FASTA filepath (e.g. 42mer_k14_probe.fa)
                                    """)
    required_arguments.add_argument('-o', '--implied_output', action='store', required=True,
                                    default=None, type=str, 
                                    help="""Specify the implied barcode output FASTA filepath (e.g. 42mer_k14_implied.fa)
                                    """)
    # required_arguments.add_argument('-c', '--combi_output', action='store', required=True,
    #                                 default=None, type=str, 
    #                                 help="""Specify the output TSV filepath with all nb-probe-implied barcode combinations (e.g. combi.tsv)
    #                                 """)
    
    args = user_input.parse_args()
    check_args(args)
    check_file_exists(args.implied_output)
    # check_file_exists(args.combi_output)
    params = {}
    params["nb_input_filepath"] = args.nb_input
    params["probe_input_filepath"] = args.probe_input
    params["implied_output_filepath"] = args.implied_output
    # params["combi_output_filepath"] = args.combi_output
    get_implied(params)

if __name__ == '__main__':
    main()