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

Lastly, from gtk-lab Amazon S3 bucket called `saber-plish-barcodes`, click into the folder `reference_genomes` and download the folder in it called `GRCh38.p14/` (with all its contents). Put the folder in the directory called `reference_genomes` which should be in the same directory as `workflow.sh`. In other words, the directory tree should be as such:
- `workflow.sh`
- `reference_genomes`
  - `GRCh38.p14`
    - `genome_assemblies_genome_fasta.tar`
    - `GRCh38_index.1.bt2`
    - `GRCh38_index.2.bt2`
    - `GRCh38_index.3.bt2`
    - `GRCh38_index.4.bt2`
    - `GRCh38_index.rev.1.bt2`
    - `GRCh38_index.rev.2.bt2`
    - `GRCh38.fa`
    - `md5checksums.txt`
    - `README.txt`
    - `report.txt`


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


# Example run-through with workflow.sh

In this section, we will show an example of running the `workflow.sh` script. All user inputs will be **bolded** while all lines printed on the commandline will be in `codeblocks`.

**bash workflow.sh** \
`STEP 1: Generating sequence libraries with SeqWalk`      
`Length (L) of the sequences in the library? (must be an even number) :` **42**  
`Length of sequences L = 42`   
`Degree of orthogonality (k) ? : ` **13**     
`Set SeqWalk RCfree argument to True or False? [t/f] (try to set 't' only for low values of L and/or k as this is computationally expensive) : ` **f**     
<br/>
`Starting SeqWalk for the following parameters:
        L=42, 
        k=13, 
        bases=ATCG, 
        prevented_patterns=['AAAA', 'TTTT', 'CCCC', 'GGGG'], 
        RCfree=False`  
<br/>
`SeqWalk completed!`  
`Size of library: 1704096`  
`Time taken: 108.13844680786133 seconds`   
<br/>
`STEP 2: Filtering base composition of sequences`  
`Lower bound for NANOBODY barcode and PROBE barcode nucleotide composition in percentage (%) (e.g. 23.5): ` **23.5**  
`Upper bound for NANOBODY barcode and PROBE barcode nucleotide composition in percentage (%) (e.g. 26.5): ` **26.5**    
`Base content filtering completed!`  
`Size of library before filtering nucleotide composition: 1704096`  
`Size of library after filtering nucleotide composition: 15787`
<br/>
`STEP 3: Align sequences to the human genome with Bowtie2 and remove any which align`  
`Is the reference genome indexes at this location: ./reference_genomes/GRCh38.p14/GRCh38_index? [y/n] :` **y**  
`Time loading reference: 00:00:00`  
`Time loading forward index: 00:00:02`  
`Time loading mirror index: 00:00:00`  
`Multiseed full-index search: 00:00:01`  
`15787 reads; of these:`  
&nbsp;&nbsp;&nbsp;&nbsp;`15787 (100.00%) were unpaired; of these:`  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;`15781 (99.96%) aligned 0 times`  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;`4 (0.03%) aligned exactly 1 time`  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;`2 (0.01%) aligned >1 times`  
`0.04% overall alignment rate`  
`Time searching: 00:00:03`  
`Overall time: 00:00:03`  
<br/>
`Number of sequences before bowtie2 alignment step: 15792`  
`Number of sequences after bowtie2 alignment step: 15781`  
<br/>
`STEP 4: Remove subsequence reverse complements from existing library`  
`Reversed complements have been removed!`  
`Size of library before RC removal: 15781`  
`Size of library after RC removal: 13085`  
<br/>
`STEP 5: Assign sequences to either NANOBODY barcode or PROBE barcode libraries (total assigned cannot exceed total number of available sequences)`  
`Number of sequences to assign as NANOBODY barcodes: ` **85**  
`Number of sequences to assign as PROBE barcodes: ` **13000**  
`From library of size 13085, assigning top 13085`  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;`sequences with highest complexity to nanobody barcode and ISH probe barcode libraries.`  
`85 sequences has been assigned to nanobody barcode libraries.`  
`13000 sequences has been assigned to ISH probe barcode libraries.`  
<br/>
`STEP 6: Generate starting IMPLIED barcode library by cross-joining NANOBODY barcode and PROBE barcode libraries`  
`Number of Implied barcodes derived: 1105000`  
<br/>
`STEP 7: Align IMPLIED barcode sequences to the human genome with Bowtie2 and remove corresponding NANOBODY and PROBE barcodes whose IMPLIED barcode aligns`  
`Time loading reference: 00:00:00`  
`Time loading forward index: 00:00:01`  
`Time loading mirror index: 00:00:01`  
`Multiseed full-index search: 00:00:21`  
`1105000 reads; of these:`  
&nbsp;&nbsp;&nbsp;&nbsp;`1105000 (100.00%) were unpaired; of these:`  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;`1104886 (99.99%) aligned 0 times`  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;`110 (0.01%) aligned exactly 1 time`  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;`4 (0.00%) aligned >1 times`  
`0.01% overall alignment rate`  
`Time searching: 00:00:23`  
`Overall time: 00:00:23`  
<br/>
`Final NANOBODY barcode library count: 69`  
`Final ISH PROBE barcode library count: 12889`  
<br/>
`STEP 8: Re-generate IMPLIED barcode library by cross-joining final NANOBODY barcode and final PROBE barcode libraries`  
`Number of Implied barcodes derived: 889341`  
<br/>
`Final NANOBODY barcode library file: data/42mer_k13/42mer_k13_nb.fa`  
`Final PROBE barcode library file: data/42mer_k13/42mer_k13_probe.fa`  
`Final IMPLIED barcode library file: data/42mer_k13/42mer_k13_implied.fa`  
`Final .tsv file with all NANOBODY-PROBE-IMPLIED barcode combinations: data/42mer_k13/42mer_k13_combi.tsv`  
