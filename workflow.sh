#!/bin/bash

printf "STEP 1: Generating sequence libraries with SeqWalk \n"
read -p "Length (L) of the sequences in the library? (must be an even number) : " L
if [[ $L%2 -eq 0 ]]; then
    printf "Length of sequences L = $L \n"
else
    echo "ERROR: L must be an even number"
    exit 125
fi

read -p "Degree of orthogonality (k) ? : " k

basename="${L}mer_k${k}"
basedir="data/${basename}"

printf "Degree of orthogonality k = $k \n"
read -p "Set SeqWalk RCfree argument to True or False? [t/f] (try to set 't' only for low values of L and/or k as this is computationally expensive) : " rcfree
if [[ $rcfree == "t" ]]; then
    if [[ $k%2 -eq 0 ]]; then
        echo "ERROR: For SeqWalk's RCfree to be True, k must be an odd number"
        exit 125
    else
        python3 generate_seqs.py -L ${L} -k ${k} -o ${basedir}/${basename}.fa -r
    fi
elif [[ $rcfree == "f" ]]; then
    python3 generate_seqs.py -L ${L} -k ${k} -o ${basedir}/${basename}.fa
else
    printf "Answer can only be 't' or 'f' (case sensitive). Answer ambiguous, defaulting to 'f'."
    rcfree="f"
    python3 generate_seqs.py -L ${L} -k ${k} -o ${basedir}/${basename}.fa
fi
printf "\n"

printf "STEP 2: Filtering base composition of sequences \n"
read -p "Lower bound for NANOBODY barcode and PROBE barcode nucleotide composition in percentage (%) (e.g. 23.5): " lower_com
read -p "Upper bound for NANOBODY barcode and PROBE barcode nucleotide composition in percentage (%) (e.g. 26.5): " upper_com
python3 filter_composition.py -f ${basedir}/${basename}.fa -o ${basedir}/${basename}_com.fa -a ${lower_com} ${upper_com} -t ${lower_com} ${upper_com} -c ${lower_com} ${upper_com} -g ${lower_com} ${upper_com}
printf "\n"

printf "STEP 3: Align sequences to the human genome with Bowtie2 and remove any which align \n"
ref_genome="./reference_genomes/GRCh38.p14/GRCh38_index"
read -p "Is the reference genome indexes at this location: ${ref_genome}? [y/n] : " ref_genome_ans
if [[ $ref_genome_ans == "n" ]]; then 
    read -p "Location of reference genome indexes for Bowtie2 (e.g. ./reference_genomes/GRCh38.p14/GRCh38_index): " ref_genome
fi
bowtie2 -x ${ref_genome} -f ${basedir}/${basename}_com.fa --no-hd -t -k 100 --very-sensitive-local -S ${basedir}/${basename}_com_bt.sam
printf "\n"
python3 get_unaligned.py -f ${basedir}/${basename}_com_bt.sam -o ${basedir}/${basename}_com_bt.fa
printf "\n"

printf "STEP 4: Remove subsequence reverse complements from existing library \n"
if [[ $rcfree == "t" ]]; then
    printf "Skipping this step as RCfree was already set to True when generating sequences with SeqWalk. \n"
elif [[ $rcfree == "f" ]]; then
    python3 remove_rc.py -f ${basedir}/${basename}_com_bt.fa -o ${basedir}/${basename}_com_bt_rcfree.fa -k ${k}
fi
printf "\n"

printf "STEP 5: Assign sequences to either NANOBODY barcode or PROBE barcode libraries (total assigned cannot exceed total number of available sequences) \n"
read -p "Number of sequences to assign as NANOBODY barcodes: " num_nb_bc
read -p "Number of sequences to assign as PROBE barcodes: " num_probe_bc

python3 assign_bc.py -f ${basedir}/${basename}_com_bt_rcfree.fa -n ${basedir}/${basename}_nb_pre.fa -p ${basedir}/${basename}_probe_pre.fa -nb ${num_nb_bc} -probe ${num_probe_bc}
printf "\n"

printf "STEP 6: Generate starting IMPLIED barcode library by cross-joining NANOBODY barcode and PROBE barcode libraries \n"
python3 get_implied.py -n ${basedir}/${basename}_nb_pre.fa -p ${basedir}/${basename}_probe_pre.fa -o ${basedir}/${basename}_implied_pre.fa -t ${basedir}/${basename}_combi_pre.tsv
printf "\n"

printf "STEP 7: Align IMPLIED barcode sequences to the human genome with Bowtie2 and remove corresponding NANOBODY and PROBE barcodes whose IMPLIED barcode aligns  \n"
bowtie2 -x ${ref_genome} -f ${basedir}/${basename}_implied_pre.fa --no-hd -t -k 100 --very-sensitive-local -S ${basedir}/${basename}_implied_pre_bt.sam
printf "\n"
python3 get_corresponding_unaligned.py -i ${basedir}/${basename}_implied_pre_bt.sam -ni ${basedir}/${basename}_nb_pre.fa -pi ${basedir}/${basename}_probe_pre.fa -no ${basedir}/${basename}_nb.fa -po ${basedir}/${basename}_probe.fa
printf "\n"

printf "STEP 8: Re-generate IMPLIED barcode library by cross-joining final NANOBODY barcode and final PROBE barcode libraries \n"
python3 get_implied.py -n ${basedir}/${basename}_nb.fa -p ${basedir}/${basename}_probe.fa -o ${basedir}/${basename}_implied.fa -t ${basedir}/${basename}_combi.tsv
printf "\n"

printf "Final NANOBODY barcode library file: ${basedir}/${basename}_nb.fa \n"
printf "Final PROBE barcode library file: ${basedir}/${basename}_probe.fa \n"
printf "Final IMPLIED barcode library file: ${basedir}/${basename}_implied.fa \n"
printf "Final .tsv file with all NANOBODY-PROBE-IMPLIED barcode combinations: ${basedir}/${basename}_combi.tsv \n"