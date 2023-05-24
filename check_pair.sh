read -p "Nanobody barcode sequence: " nb_seq
read -p "Probe barcode sequence: " probe_seq

python3 get_implied_simple.py -n $nb_seq -p $probe_seq -o check/combi.fa
printf "\n"

ref_genome="./reference_genomes/GRCh38.p14/GRCh38_index"
read -p "Is the reference genome indexes at this location: ${ref_genome}? [y/n] : " ref_genome_ans
if [[ $ref_genome_ans == "n" ]]; then 
    read -p "Location of reference genome indexes for Bowtie2 (e.g. ./reference_genomes/GRCh38.p14/GRCh38_index): " ref_genome
fi

bowtie2 -x ${ref_genome} -f check/combi.fa --no-hd -t -k 100 --very-sensitive-local -S check/combi.sam
printf "Alignment completed, inspect file at check/combi.sam for alignment with genome."