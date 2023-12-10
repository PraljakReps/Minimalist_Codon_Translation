#!/usr/bin/env sh

export DIR="$(dirname "$(pwd)")"

export input_fasta_path='/project/andrewferguson/niksapraljak/Minimalist_Codon_Translation/outputs/output_SH3_DNA_test.fasta'
export output_fasta_path='/project/andrewferguson/niksapraljak/Minimalist_Codon_Translation/outputs/output_SH3_optimized_DNA_test.fasta'
export lim_GC_high=0.6 # default: 0.60
export lim_GC_low=0.45 # default: 0.45
export codon_bias_path='/project/andrewferguson/niksapraljak/Minimalist_Codon_Translation/Codon_bias/YeastProteomeContent.tsv' # codon frequency file (depends on specie)
export seed=42

python ../main.py \
	--input ${input_fasta_path} \
	--output ${output_fasta_path} \
	--lim_GC_high ${lim_GC_high} \
	--lim_GC_low ${lim_GC_low} \
	--codon_bias_path ${codon_bias_path} \
	--seed ${seed} \
