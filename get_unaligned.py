import argparse
from pathlib import Path

def check_args(args):
    if args.file == args.output:
        raise Exception("Input and output file cannot be the same.")
    if not args.file.endswith(".sam"):
        raise Exception("Input file must be a SAM file (.sam).")
    
def check_file_exists(output_filepath):
    file = Path(output_filepath)
    if file.is_file():
        reply = input(f"File {output_filepath} already exists. Overwrite? [y/n] ")
        if reply == "y":
            return
        else:
            print("Command terminated.")
            exit()

def get_unaligned(params):
    with open(params.get("input_filepath"), "r") as input_filepath:
        with open(params.get("output_filepath"), "w") as output_filepath:
            input_count = 0
            count = 0
            for line in input_filepath:
                input_count += 1
                line_list = list(line.strip("\n").split("\t"))
                id = line_list[0]
                flag = line_list[1]
                sequence = line_list[9]
                if flag == "4":
                    output_filepath.write(f">{id}\n{sequence}\n")
                    count += 1
    print(f"Number of sequences before bowtie2 alignment step: {input_count}")
    print(f"Number of sequences after bowtie2 alignment step: {count}")

def main():
    user_input = argparse.ArgumentParser()
    required_arguments = user_input.add_argument_group('required arguments')
    required_arguments.add_argument('-f', '--file', action='store', default=None, type=str, required=True,
                                    help='Path to input FASTA file')
    required_arguments.add_argument('-o', '--output', action='store', required=True,
                                    default=None, type=str, 
                                    help="""Specify the output FASTA filepath (e.g. my_sequeces_filtered.fa)
                                    """)
    
    args = user_input.parse_args()
    check_args(args)
    check_file_exists(args.output)
    params = {}
    params["input_filepath"] = args.file
    params["output_filepath"] = args.output
    get_unaligned(params)

if __name__ == '__main__':
    main()