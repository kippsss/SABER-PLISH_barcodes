# SABER-PLISH_barcodes

This repository contains all the scripts used in the computational design of DNA barcode libraries for use in SABER-PLISH.

The DNA barcode libraries generated are not uploaded here due to large file sizes. They will instead be uploaded to a private Amazon S3 bucket.

# Requirements to satisfy

Firstly, ensure that you have all the python dependencies required for the entire workflow. This is listed in `requirements.txt`. To install them, simply input `pip install -r requirements.txt` into the commandline.

Secondly, ensure that you have Bowtie2 installed. Refer to Bowtie2 manual https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml or this anaconda/bioconda page https://anaconda.org/bioconda/bowtie2 to get started on the installation.

Thirdly, ensure that you have at least the following scripts in the same directory:
- `assign_bc.py`
- `filter_composition.py`
- `generate_seqs.py`
- `get_corresponding_unaligned.py`
- `get_implied.py`
- `get_unaligned.py`
- `remove_rc.py`
- `workflow.sh`

# Using workflow.sh to generate barcode libraries

1) In the commandline, run `bash workflow.sh`.
2) A prompt `Length (L) of the sequences in the library? (must be an even number)` will be shown.
3) A prompt `Degree of orthogonality (k) ?` will be shown. It is recommended that k <= 0.4L - 0.5L.
4) A prompt `Set SeqWalk RCfree argument to True or False? [t/f] (try to set 't' only for low values of L and/or k as this is computationally expensive)` will be shown. For high values of L or k, you are recommended to set this to False (e.g. L=42 and k=13 with RCfree=True is too computationally expensive to compute).
5) SeqWalk will start generating the starting sequence library with the given input parameters.
6) A prompt `Lower bound for NANOBODY barcode and PROBE barcode nucleotide composition in percentage (%) (e.g. 23.5): ` will be shown. This is the value *a*, such that for all sequences, a type of nucleotide (A/T/C/G) will compose >= *a*% of the sequence.
7) A prompt `Upper bound for NANOBODY barcode and PROBE barcode nucleotide composition in percentage (%) (e.g. 26.5): ` will be shown. This is the value *b*, such that for all sequences, a type of nucleotide (A/T/C/G) will compose <= *b*% of the sequence.
8) A nucleotide composition filter algorithm will run, such that only all sequences that do not fall within the range *a*% <= x <= *b*% will be removed.
9) A prompt `Is the reference genome indexes at this location: ./reference_genomes/GRCh38.p14/GRCh38_index? [y/n] : ` will be shown. It is ideal that you have downloaded the GRCh38 index files found in the Amazon S3 bucket, into the location stated in the prompt. If so, simply input `y`, else, input `n`. If `n` is received, another prompt `Location of reference genome indexes for Bowtie2 (e.g. ./reference_genomes/GRCh38.p14/GRCh38_index): ` will be shown. Input accordingly.
10) Bowtie2 will be used to align sequences. Any sequence with alignment to the genome will be removed.
11) A prompt `Number of sequences to assign as NANOBODY barcodes: ` will be shown. Input desired number *n* of nanobody barcodes.
12) A prompt `Number of sequences to assign as PROBE barcodes: ` will be shown. Input desired number *p* of probe barcodes. (note that *n* + *p* cannot exceed number of remainding sequences after Bowtie2 filtering.
13) The starting implied barcode library will generated from a cross-join of the assigned nanobody barcode library and probe barcode library (size of implied barcode library = *n* x *p*)
14) Starting implied barcode library will be aligned to the genome with Bowtie2. For every aligned implied barcode, its corresponding nanobody and probe barcodes will be removed from their libraries.
15) The new nanobody and probe barcode libraries (which are now final) are re-cross-joined to form the new and final implied barcode library.
16) Combinations of all nanobody-probe-implied trios will be generated as a .tsv file.
17) The file names of the final libraries will be printed. (e.g. `Final IMPLIED barcode library file: data/42mer_k12/42mer_k12_implied.fa`)
18) The file name of the nanobody-probe-implied trio combinations will also be printed (e.g. `Final .tsv file with all NANOBODY-PROBE-IMPLIED barcode combinations: data/42mer_k12/42mer_k12_combi.tsv`)
