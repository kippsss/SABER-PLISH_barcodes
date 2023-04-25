from seqwalk import design
import argparse
from pathlib import Path
import os
import time

def get_output_filepath(output_filepath, L, k):
    if output_filepath is None:
        output_filename = f"{L}mer_k{k}.fa"
        output_filepath = "./" + output_filename
    else:
        pass
    return output_filepath

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

def write_comments(params, output_filepath, exclude=list()):
    for param_key, param_value in params.items():
        if param_key in exclude:
            continue
        else:
            output_filepath.write(f";{param_key}={param_value}\n")
    output_filepath.write("\n")

def initiate_seqwalk(params):

    start = time.time()
    
    print(
    f"""
    Starting SeqWalk for the following parameters:
        L={params.get("L")}, 
        k={params.get("k")}, 
        bases={params.get("bases")}, 
        prevented_patterns={params.get("prevented_patterns")}, 
        RCfree={params.get("rcfree")}
    """
    )

    library_max_size = design.max_size(params.get("L"), 
                                    params.get("k"), 
                                    alphabet=params.get("bases"), 
                                    prevented_patterns=params.get("prevented_patterns"), 
                                    RCfree=params.get("rcfree"))

    with open(params.get("output_filepath"), "w") as output_filepath:
        write_comments(params, output_filepath, exclude=["output_filepath"])
        count = 0
        L = params.get("L")
        k = params.get("k")
        for seq in library_max_size:
            count += 1
            output_filepath.write(f">seq{count}\n{seq}\n")

    end = time.time()
    total_time = end - start

    print(f"""SeqWalk completed!\nSize of library: {count}""")
    print(f"Time taken: {total_time} seconds")
    

def main():
    
    user_input = argparse.ArgumentParser()
    required_arguments = user_input.add_argument_group('required arguments')
    optional_arguments = user_input.add_argument_group('optional arguments')
    required_arguments.add_argument('-L', '--L', action='store', default=None, type=int, required=True,
                                    help='L refers to the desired length of the sequences.')
    required_arguments.add_argument('-k', '--k', action='store', default=None, type=int, required=True,
                                    help="""k refers to the orthoganality of the sequence library where
                                    there will be no subsequence of length >=k that appears more than once in
                                    the library of sequences.""")
    optional_arguments.add_argument('-b', '--bases', action='store', 
                                    default="ATCG", type=str, 
                                    help="""bases is a string where each character is a base allowed for use in
                                    generating the sequence library.
                                    Default="ATCG"
                                    """)
    optional_arguments.add_argument('-p', '--prevented_patterns', action='store', nargs='+', 
                                    default=["AAAA","TTTT","CCCC","GGGG"], type=str, 
                                    help="""Please specify subsstrings delimited by a space where each string is a pattern
                                    that must be prevented from occuring in the library. 
                                    Default=AAAA TTTT CCCC GGGG""")
    optional_arguments.add_argument('-r', '--rcfree', action='store_true', 
                                    help="""If -r is specified, there will be no reverse complement for any subsequence
                                    of length k, in that case, k must be an odd number for SeqWalk to work.
                                    Not recommended if k >= 13 due to high time complexity.
                                    Default=False""")
    optional_arguments.add_argument('-o', '--output', action='store', 
                                    default=None, type=str, 
                                    help="""Specify the output filepath (e.g. my_sequences.fa)
                                    Default="\{L\}mer_k\{k\}.fa in current directory where generate_seqs.py is run."
                                    """)
    args = user_input.parse_args()
    output_filepath = get_output_filepath(args.output, args.L, args.k)
    check_file_exists(output_filepath)
    params = {
        "L": args.L,
        "k": args.k,
        "bases": args.bases,
        "prevented_patterns": args.prevented_patterns,
        "rcfree": args.rcfree,
        "output_filepath": output_filepath,
    }
    initiate_seqwalk(params)

if __name__ == '__main__':
    main()
