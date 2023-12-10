


import argparse
import source.DNA_tools as dt
import source.robustness_tools_072921 as rt # what is this? 
import random
import sys
import source.DNA_assembly_tools as assembly_dt
import source.reverse_translation_tools as reverse_tt



def set_seed(args: any) -> None:
    random.seed(args.seed)
    return None


def get_args(parser: any) -> None:
    
    parser.add_argument('--input', default=None, type=str,
                            help='fasta with your protein sequence of interest.')
    parser.add_argument('--output_dna_path', default=None, type=str,
                            help='path with your dna sequence.')
    parser.add_argument('--codon_bias_path', default=None, type=str,
                            help='path to the codon frequency.')
    parser.add_argument('--seed', default=42, type=int,                                                                                   
                            help='random seed for reproducibility')
    parser.add_argument('--aa_header', default='amino_acid', type=str,                                                                                   
                            help='amino acid column name')
    parser.add_argument('--codon_header', default='codon', type=str,                                                                                   
                            help='codon column name')
    parser.add_argument('--frequency_header', default='norm_freq', type=str,                                                                                   
                            help='codon frequency column name')
    
    return parser



#################
# Main function #
#################

def run_reverse_translation(args: any) -> None:

    # read fasta or csv spreadsheet
    protein_seqs = reverse_tt.read_proteins(args=args)
    print(protein_seqs)

    # convert list of aa sequences to nt sequences
    dna_seqs = reverse_tt.convert_aa_2_nt(
            args=args,
            protein_sequences=protein_seqs
    )
    print(dna_seqs)

    # save fasta, csv, or both.
    reverse_tt.save_dna_sequences(
            args=args,
            dna_sequences=dna_seqs
    )


    return 
    

if __name__ == '__main__':
    
    
    # get variables
    parser = argparse.ArgumentParser(description='AA Reverse Translation and Codon Optimization')
    parser = get_args(parser=parser)
    args = parser.parse_args()

    # for reproducibility
    set_seed(args=args)
    
    # run DNA assembly and Codon optimizatipm
    run_reverse_translation(args=args)


    

