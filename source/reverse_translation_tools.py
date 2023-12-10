from Bio import SeqIO
import csv
import random


def read_proteins(args: any):
    # return protein design library as a dictionary
    return {record.id: str(record.seq) for record in SeqIO.parse(args.input, 'fasta')}


def parse_codon_usage(args: any):
    
    # codon spreadsheet path
    filename = args.codon_bias_path

    codon_usage = {}
    # Determine the delimiter based on the file extension
    delimiter = '\t' if filename.endswith('.tsv') else ','

    with open(filename, 'r') as csvfile:

        reader = csv.DictReader(csvfile, delimiter=delimiter)

        for row in reader:

            amino_acid = row[args.aa_header]
            codon = row[args.codon_header]
            frequency = float(row[args.frequency_header])

            if amino_acid not in codon_usage:
                codon_usage[amino_acid] = []

            codon_usage[amino_acid].append((codon, frequency))

    return codon_usage


def weighted_choice(choices):

    total = sum(weight for codon, weight in choices)

    r = random.uniform(0, total)

    upto = 0
    for codon, weight in choices:

        if upto + weight > r:
            return codon
        upto += weight


def reverse_translate(protein_sequence, codon_usage):

    dna_sequence = ''
    for amino_acid in protein_sequence:

        codons = codon_usage.get(amino_acid, [])
        if not codons:
            raise ValueError(f"No codons found for amino acid {amino_acid}")
        
        codon = weighted_choice(codons)
        dna_sequence += codon
    return dna_sequence


def convert_aa_2_nt(
        args: any,
        protein_sequences: dict
    ) -> dict:

    # get codon frequencies 
    codon_usage = parse_codon_usage(args=args)
    
    # prepare an empty dna dictionary 
    dna_sequences = {}

    for  seq_name, seq in protein_sequences.items():

        dna_sequence = reverse_translate(
                protein_sequence=seq,
                codon_usage=codon_usage
        )

        dna_sequences[seq_name] = dna_sequence


    return dna_sequences




def save_dna_sequences(
        args: any,
        dna_sequences: dict
    ) -> None:

    # Determine the file extension
    filename = args.output_dna_path
    file_extension = filename.split('.')[-1]

    if file_extension == 'fasta':
        with open(filename, 'w') as file:
            for seq_name, dna_seq in dna_sequences.items():
                file.write(f">{seq_name}\n{dna_seq}\n")
    elif file_extension in ['csv', 'tsv']:
        delimiter = ',' if file_extension == 'csv' else '\t'
        with open(filename, 'w') as file:
            for seq_name, dna_seq in dna_sequences.items():
                file.write(f"{seq_name}{delimiter}{dna_seq}\n")


