#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 18 13:16:01 2021

@author: bryanandrews1
@editor: email -> niksapraljak1@gmail.com | Github -> PraljakReps
"""

from optparse import OptionParser
import DNA_tools as dt
import robustness_tools_072921 as rt # what is this? 
import random
import sys

def set_seed(args: any) -> None:
    random.seed(args.seed)
    return None

def get_args(parser: any) -> None:
    
    parser.add_argument('--input', default=None, type=str,
            help='fasta with your protein sequence of interest')
    parser.add_argument('--output', default=None, type=str,
            help='fasta with your optimized sequence')

    return parser


#################
# Main function #
#################

def run_main(args: any) -> None:

    # get protein sequence and DNA sequences
    seq_name, DNA_seq = dt.read_fasta(args.input)
    



    


if __name__ == "__main__":


    # get variables
    parser = argparse.ArgumentParser(description='AA Reverse Translation and Codon Optimization')
    parser = get_args(parser=parser)
    args = parser.parse_args()

    # for reproducibility
    set_seed(args=args)

        
    # run DNA assembly and Codon optimizatipm
    run_main(args=args)






