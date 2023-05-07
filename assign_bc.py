import argparse
from pathlib import Path
import math
import random
from Bio.Seq import Seq
from Bio import SeqIO

def check_args(args):
    if args.file == args.nb_output or args.file == args.probe_output:
        raise Exception("Input and output files cannot be the same.")
    if not args.file.endswith(".fa") and not args.file.endswith(".fasta"):
        raise Exception("Input file must be a FASTA file (.fa).")
    if not args.nb_output.endswith(".fa") and not args.nb_output.endswith(".fasta"):
        raise Exception("Nanobody barcodes output file must be a FASTA file (.fa).")
    if not args.probe_output.endswith(".fa") and not args.probe_output.endswith(".fasta"):
        raise Exception("ISH probe barcodes output file must be a FASTA file (.fa).")
    
def check_file_exists(output_filepath):
    file = Path(output_filepath)
    if file.is_file():
        reply = input(f"File {output_filepath} already exists. Overwrite? [y/n] ")
        if reply == "y":
            return
        else:
            print("Command terminated.")
            exit()

def get_complexity_numerator(seq):
    length = len(seq)
    num_words = []
    for word_size in range(1, length):
        words = {}
        for i in range(length-word_size+1):
            word = seq[i: i+word_size]
            if word not in words.keys():
                words[word] = 1
        num_words.append(len(words))
    complexity_numerator = math.prod(num_words)
    return complexity_numerator

def check_bc_availability(input_count, num_bc_required):
    if num_bc_required > input_count:
        raise Exception(f"""Number of barcodes requested is more than number of barcodes available.
        \nNumber of barcodes requested: {num_bc_required}\nNumber of barcodes available: {input_count}""")

def assign_bc(params):
    input_fasta = SeqIO.parse(open(params.get("input_filepath")),'fasta')
    input_count = 0
    nb_bc_count = 0
    probe_bc_count = 0
    complexity_dict = {}
    for line in input_fasta:
        input_count += 1
        id, seq = line.id, str(line.seq)
        complexity_numerator = get_complexity_numerator(seq)
        complexity_dict[(id, seq)] = complexity_numerator
    sorted_complexity_dict = dict(sorted(complexity_dict.items(), key=lambda x:x[1], reverse=True))

    num_nb_bc = params.get("num_nb_barcodes")
    num_probe_bc = params.get("num_probe_barcodes")
    num_bc_required = num_nb_bc + num_probe_bc
    check_bc_availability(input_count, num_bc_required)

    top_complex_bc = list(sorted_complexity_dict.keys())[:num_bc_required]
    random.shuffle(top_complex_bc)
    
    # Take first n sequences for nanobody barcodes, and remainding sequences after n for probe barcodes
    nb_bc_list = top_complex_bc[:num_nb_bc]
    probe_bc_list = top_complex_bc[num_nb_bc:]

    print(f"""From library of size {input_count}, assigning top {num_bc_required} 
    sequences with highest complexity to nanobody barcode and ISH probe barcode libraries.""")

    with open(params.get("nb_output_filepath"), "w") as nb_output_filepath:
        for nb_bc in nb_bc_list:
            nb_bc_id, nb_bc_seq = nb_bc
            nb_output_filepath.write(f">nb-{nb_bc_id}\n{nb_bc_seq}\n")
            nb_bc_count += 1
    print(f"{nb_bc_count} sequences has been assigned to nanobody barcode libraries.")

    with open(params.get("probe_output_filepath"), "w") as probe_output_filepath:
        for probe_bc in probe_bc_list:
            probe_bc_id, probe_bc_seq = probe_bc
            probe_output_filepath.write(f">probe-{probe_bc_id}\n{probe_bc_seq}\n")
            probe_bc_count += 1
    print(f"{probe_bc_count} sequences has been assigned to ISH probe barcode libraries.")

def main():
    user_input = argparse.ArgumentParser()
    required_arguments = user_input.add_argument_group('required arguments')
    required_arguments.add_argument('-f', '--file', action='store', default=None, type=str, required=True,
                                    help='Path to input FASTA file')
    required_arguments.add_argument('-n', '--nb_output', action='store', required=True,
                                    default=None, type=str, 
                                    help="""Specify the output FASTA filepath for NANOBODY barcodes(e.g. 42mer_k14_nb.fa)
                                    """)
    required_arguments.add_argument('-p', '--probe_output', action='store', required=True,
                                    default=None, type=str, 
                                    help="""Specify the output FASTA filepath for PROBE barcodes (e.g. 42mer_k14_probe.fa)
                                    """)
    required_arguments.add_argument('-nb', '--num_nb_barcodes', action='store', required=True,
                                    default=100, type=int, 
                                    help="""Specify number of nanobody barcodes (e.g. 100)
                                    """)
    required_arguments.add_argument('-probe', '--num_probe_barcodes', action='store', required=True,
                                    default=20000, type=int, 
                                    help="""Specify number of probe_barcodes barcodes (e.g. 20000)
                                    """)
    
    args = user_input.parse_args()
    check_args(args)
    check_file_exists(args.nb_output)
    check_file_exists(args.probe_output)
    params = {}
    params["input_filepath"] = args.file
    params["nb_output_filepath"] = args.nb_output
    params["probe_output_filepath"] = args.probe_output
    params["num_nb_barcodes"] = args.num_nb_barcodes
    params["num_probe_barcodes"] = args.num_probe_barcodes
    assign_bc(params)

if __name__ == '__main__':
    main()