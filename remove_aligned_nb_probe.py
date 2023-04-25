import argparse
from pathlib import Path
import math
import random
from Bio.Seq import Seq
from Bio import SeqIO

def check_args(args):
    # if args.nb_input == args.implied_output or args.probe_input == args.implied_output:
    #     raise Exception("Input and output files cannot be the same.")
    if not args.implied_bt_input.endswith(".sam"):
        raise Exception("Aligned implied barcode input file must be SAM format (.sam)")
    if not args.nb_input.endswith(".fa") and not args.nb_input.endswith(".fasta"):
        raise Exception("Nanobody barcode input file must be a FASTA file (.fa).")
    if not args.probe_input.endswith(".fa") and not args.probe_input.endswith(".fasta"):
        raise Exception("Probe barcode input file must be a FASTA file (.fa).")
    if not args.nb_output.endswith(".fa") and not args.nb_output.endswith(".fasta"):
        raise Exception("Nanobody barcode output file must be a FASTA file (.fa).")
    if not args.probe_output.endswith(".fa") and not args.probe_output.endswith(".fasta"):
        raise Exception("ISH probe barcode output file must be a FASTA file (.fa).")

def check_file_exists(output_filepath):
    file = Path(output_filepath)
    if file.is_file():
        reply = input(f"File {output_filepath} already exists. Overwrite? [y/n] ")
        if reply == "y":
            return
        else:
            print("Command terminated.")
            exit()

def remove_aligned_nb_probe(params):
    with open(params.get("implied_bt_input_filepath"), "r") as implied_bt_input_filepath:
        nb_to_remove = []
        probe_to_remove = []
        for line in implied_bt_input_filepath:
            line_list = list(line.strip("\n").split("\t"))
            implied_id = line_list[0]
            flag = line_list[1]
            if flag != "4":
                nb_id, probe_id = implied_id.split("_")
                if nb_id not in nb_to_remove:
                    nb_to_remove.append(nb_id)
                if probe_id not in probe_to_remove:
                    probe_to_remove.append(probe_id)

    with open(params.get("nb_output_filepath"), "w") as nb_output_filepath:
        final_nb_count = 0
        nb_input_fasta = SeqIO.parse(open(params.get("nb_input_filepath")),'fasta')
        for nb_line in nb_input_fasta:
            nb_id, nb_seq = nb_line.id, str(nb_line.seq)
            if nb_id not in nb_to_remove:
                nb_output_filepath.write(f">{nb_id}\n{nb_seq}\n")
                final_nb_count += 1

    with open(params.get("probe_output_filepath"), "w") as probe_output_filepath:
        final_probe_count = 0
        probe_input_fasta = SeqIO.parse(open(params.get("probe_input_filepath")),'fasta')
        for probe_line in probe_input_fasta:
            probe_id, probe_seq = probe_line.id, str(probe_line.seq)
            if probe_id not in probe_to_remove:
                probe_output_filepath.write(f">{probe_id}\n{probe_seq}\n")
                final_probe_count += 1

    print(f"Final NANOBODY barcode pool count: {final_nb_count}")
    print(f"Final ISH PROBE barcode pool count: {final_probe_count}")

def main():
    user_input = argparse.ArgumentParser()
    required_arguments = user_input.add_argument_group('required arguments')
    required_arguments.add_argument('-i', '--implied_bt_input', action='store', required=True,
                                    default=None, type=str, 
                                    help="""Specify the SAM file after alignment of implied barcodes to genome with Bowtie2 (42mer_k14_implied_bt.sam)
                                    """)
    required_arguments.add_argument('-ni', '--nb_input', action='store', required=True,
                                    default=None, type=str, 
                                    help="""Specify the nanobody barcode input FASTA filepath (e.g. 42mer_k14_nb.fa)
                                    """)
    required_arguments.add_argument('-pi', '--probe_input', action='store', required=True,
                                    default=None, type=str, 
                                    help="""Specify the ISH probe barcode input FASTA filepath (e.g. 42mer_k14_probe.fa)
                                    """)
    required_arguments.add_argument('-no', '--nb_output', action='store', required=True,
                                    default=None, type=str, 
                                    help="""Specify the nanobody barcode output FASTA filepath after removal (e.g. 42mer_k14_nb_final.fa)
                                    """)
    required_arguments.add_argument('-po', '--probe_output', action='store', required=True,
                                    default=None, type=str, 
                                    help="""Specify the ISH probe barcode output FASTA filepath after removal (e.g. 42mer_k14_probe_final.tsv)
                                    """)
    
    args = user_input.parse_args()
    check_args(args)
    check_file_exists(args.nb_output)
    check_file_exists(args.probe_output)
    params = {}
    params["implied_bt_input_filepath"] = args.implied_bt_input
    params["nb_input_filepath"] = args.nb_input
    params["probe_input_filepath"] = args.probe_input
    params["nb_output_filepath"] = args.nb_output
    params["probe_output_filepath"] = args.probe_output
    remove_aligned_nb_probe(params)

if __name__ == '__main__':
    main()