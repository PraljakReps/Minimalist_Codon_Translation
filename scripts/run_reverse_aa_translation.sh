#!/usr/bin/env sh

export DIR="$(dirname "$(pwd)")"

export input_fasta_path='/project/andrewferguson/niksapraljak/Minimalist_Codon_Translation/inputs/SH3_test.fasta'
export output_fasta_path='/project/andrewferguson/niksapraljak/Minimalist_Codon_Translation/outputs/output_SH3_DNA_test.fasta' # take 'fasta', 'csv', or 'tsv

export codon_bias_path='/project/andrewferguson/niksapraljak/Minimalist_Codon_Translation/Codon_bias/YeastProteomeContent.tsv' # codon frequency file (depends on specie)
export seed=42
export aa_header='amino_acid'
export codon_header='codon'
export frequency_header='norm_freq'


python ../reverse_translation_main.py \
	--input ${input_fasta_path} \
	--output_dna_path ${output_fasta_path} \
	--codon_bias_path ${codon_bias_path} \
	--seed ${seed} \
	--aa_header ${aa_header} \
	--codon_header ${codon_header} \
	--frequency_header ${frequency_header} \
