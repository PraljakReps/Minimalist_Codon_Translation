#!/usr/bin/env sh

export DIR="$(dirname "$(pwd)")"

export input_fasta_path=''
export output_fasta_path=''
export lim_GC_high=0.6 # default: 0.60
export lim_GC_low=0.45 # default: 0.45
export codon_bias_path='' # codon frequency file (depends on specie)


python ../optimize_DNA_for_assembly.py \
	--input_fasta_path ${input_fasta_path} \
	--output_fasta_path ${output_fasta_path} \
	--lim_GC_high ${lim_GC_high} \
	--lim_GC_low ${lim_GC_low} \
	--codon_bias_path ${codon_bias_path} \

