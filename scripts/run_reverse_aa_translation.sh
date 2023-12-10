#!/usr/bin/env sh

export DIR="$(dirname "$(pwd)")" # Parent directory of the current working directory
export input_fasta_path='/project/andrewferguson/niksapraljak/Minimalist_Codon_Translation/inputs/SH3_test.fasta' # Path to input FASTA file
export output_fasta_path='/project/andrewferguson/niksapraljak/Minimalist_Codon_Translation/outputs/output_SH3_DNA_test.fasta' # Path to output file (fasta, csv, tsv)
export codon_bias_path='/project/andrewferguson/niksapraljak/Minimalist_Codon_Translation/Codon_bias/YeastProteomeContent.tsv' # Path to codon bias file
export seed=42 # Seed value for reproducibility in random processes
export aa_header='amino_acid' # Header for amino acid in output
export codon_header='codon' # Header for codon in output
export frequency_header='norm_freq' # Header for frequency in output


python ../reverse_translation_main.py \
	--input ${input_fasta_path} \
	--output_dna_path ${output_fasta_path} \
	--codon_bias_path ${codon_bias_path} \
	--seed ${seed} \
	--aa_header ${aa_header} \
	--codon_header ${codon_header} \
	--frequency_header ${frequency_header} \
