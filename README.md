# SABER-PLISH_barcodes

This repository contains all the scripts used in the computational design of DNA barcode libraries for use in SABER-PLISH.

Documentation to improve its reusability has yet to be completed.

The DNA barcode libraries generated are not uploaded here due to large file sizes. They will instead be uploaded to a private Amazon S3 bucket.

# Using workflow.sh to generate (1) nanobody, (2) probe and (3) implied barcode libraries

1) In the commandline, run `bash workflow.sh`.
2) A prompt `Length (L) of the sequences in the library? (must be an even number)` will be shown.
3) A prompt `Degree of orthogonality (k) ?` will be shown. It is recommended that k <= 0.4L - 0.5L.
4) A prompt `Set SeqWalk RCfree argument to True or False? [t/f] (try to set 't' only for low values of L and/or k as this is computationally expensive)` will be shown. For high values of L or k, you are recommended to set this to False (e.g. L=42 and k=13 with RCfree=True is too computationally expensive to compute).
5) SeqWalk will start generating the starting sequence library with the given input parameters.
6) A prompt `Lower bound for NANOBODY barcode and PROBE barcode nucleotide composition in percentage (%) (e.g. 23.5): ` will be shown. This is the value *a*, such that for all sequences, a type of nucleotide (A/T/C/G) will compose >= *a*% of the sequence.
7) A prompt `Upper bound for NANOBODY barcode and PROBE barcode nucleotide composition in percentage (%) (e.g. 26.5): ` will be shown. This is the value *b*, such that for all sequences, a type of nucleotide (A/T/C/G) will compose <= *b*% of the sequence.
8) A nucleotide composition filter algorithm will run, such that only all sequences that do not fall within the range *a*% <= x <= *b*% will be removed.
9) A prompt `Is the reference genome indexes at this location: ./reference_genomes/GRCh38.p14/GRCh38_index? [y/n] : ` will be shown. It is ideal that you have downloaded the GRCh38 index files found in the Amazon S3 bucket, into the location stated in the prompt. If so, simply input `y`, else, input `n`. If `n` is received, another prompt `Location of reference genome indexes for Bowtie2 (e.g. ./reference_genomes/GRCh38.p14/GRCh38_index): ` will be shown. Input accordingly.
