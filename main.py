#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 18 13:16:01 2021

@author: bryanandrews1
@editor: email -> niksapraljak1@gmail.com | Github -> PraljakReps
"""

from optparse import OptionParser
import argparse
import source.DNA_tools as dt
import source.robustness_tools_072921 as rt # what is this? 
import random
import sys
import source.DNA_assembly_tools as assembly_dt
from tqdm import tqdm


def set_seed(args: any) -> None:
    random.seed(args.seed)
    return None

def get_args(parser: any) -> None:
    
    parser.add_argument('--input', default=None, type=str,
            help='fasta with your protein sequence of interest.')
    parser.add_argument('--output', default=None, type=str,
            help='fasta with your optimized sequence.')
    parser.add_argument('--codon_bias_path', default=None, type=str,
            help='path to the codon frequency.')
    parser.add_argument('--lim_GC_high', default=0.6, type=float,
            help='GC content too high limit.')
    parser.add_argument('--lim_GC_low', default=0.45, type=float,
            help='GC content too low limit.')
    parser.add_argument('--seed', default=42, type=int,
            help='random seed for reproducibility')


    return parser

#################
# Main function #
#################

def run_main(args: any) -> None:

    # get protein sequence and DNA sequences...
    # seq_name, DNA_seq = dt.read_fasta(args.input)
    DNA_seq_dict = dt.read_fasta(fasta_file=args.input)
    
    # create an empty dict to store the final optimized DNA library
    final_DNA_seq_dict = {} 

    print("\nProcess all sequences within the input file...")
    for seq_name, DNA_seq in tqdm(DNA_seq_dict.items()):

        
        #print('Pre-optimization:\n', seq_name, DNA_seq)

        # optimize Codons based on specie frequency...
        DNA_seq = assembly_dt.codon_optimize(
                DNA_seq=DNA_seq,
                codon_bias=args.codon_bias_path
        )

        # adjust the G-C content of a given DNA sequence...
        # ensure it falls within the specified limits lim_GC_high and lim_GC_low
        DNA_seq = assembly_dt.fix_GC(
                DNA_seq=DNA_seq,
                lim_high=args.lim_GC_high,
                lim_low=args.lim_GC_low
        )

        # process a DNA sequence to disrupt microhomologoues sequences
        # microhomology -> refers to the occurence of short, repeated DNA sequences within a longer DNA sequence
        DNA_seq = assembly_dt.strip_microhomology(DNA_seq=DNA_seq)

        # process DNA sequence and modify it by identifying and altering mononucleotide tracts...
        # mononucleoditde tract -> sequence where the same nucleotide is repeated multiple times in a row (e.g. "AAAAA", "CCCCC", etc.)...
        DNA_seq = assembly_dt.strip_mononucleotide_tracts(DNA_seq=DNA_seq)

        # process DNA sequence and modify it in a way that breaks up any occurences of the triplet "GGG" (i.e. remove G-quadruplex)...
        DNA_seq = assembly_dt.strip_G_quadruplexes(DNA_seq=DNA_seq)
        
        final_DNA_seq_dict[seq_name] = DNA_seq

        #print('Post-optimization:\n', seq_name, DNA_seq)
    with open(args.output, 'w') as fasta_out:
        
        for seq_name, DNA_seq in DNA_seq_dict.items():
            fasta_out.write(">%s\n%s\n" % (seq_name, DNA_seq))

    print('Finished running AA reverse translation, DNA assembly, and Codon optimization...')
    


if __name__ == "__main__":


    # get variables
    parser = argparse.ArgumentParser(description='AA Reverse Translation and Codon Optimization')
    parser = get_args(parser=parser)
    args = parser.parse_args()

    # for reproducibility
    set_seed(args=args)
        
    # run DNA assembly and Codon optimizatipm
    run_main(args=args)






