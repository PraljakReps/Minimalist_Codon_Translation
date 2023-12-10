#!/usr/bin/env sh

export DIR="$(dirname "$(pwd)")" # Parent directory of the current working directory
export input_fasta_path='/project/andrewferguson/niksapraljak/Minimalist_Codon_Translation/outputs/output_SH3_single_DNA_test.fasta' # Path to input FASTA file
export output_fasta_path='/project/andrewferguson/niksapraljak/Minimalist_Codon_Translation/outputs/output_SH3_single_optimized_DNA_test.fasta' # Path to output FASTA file
export lim_GC_high=0.6 # Upper limit for GC content
export lim_GC_low=0.45 # Lower limit for GC content
export codon_bias_path='/project/andrewferguson/niksapraljak/Minimalist_Codon_Translation/Codon_bias/YeastProteomeContent.tsv' # Path to codon bias file
export seed=42 # Seed value for reproducibility in random processes


python ../main.py \
	--input ${input_fasta_path} \
	--output ${output_fasta_path} \
	--lim_GC_high ${lim_GC_high} \
	--lim_GC_low ${lim_GC_low} \
	--codon_bias_path ${codon_bias_path} \
	--seed ${seed} \
