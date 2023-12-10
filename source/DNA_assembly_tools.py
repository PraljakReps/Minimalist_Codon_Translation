"""
Created on Mon Oct 18 13:16:01 2021

@author: bryanandrews1
"""

import random
import source.DNA_tools as dt 
import source.robustness_tools_072921 as rt

######################
# Codon optimization #
######################

def codon_optimize(DNA_seq, codon_bias):
    codon_freqs = rt.read_codon_freqs(codon_bias)
    pro_seq = dt.translate(DNA_seq)
    if pro_seq.endswith("*"):
        syn_DNA_seq = rt.reverse_translate(pro_seq.strip("*"), CB = codon_freqs) + "TAA"
    else:
        syn_DNA_seq = rt.reverse_translate(pro_seq, CB = codon_freqs)
    return(syn_DNA_seq)


################
# Remove "GGG" #
################

def strip_G_quadruplexes(DNA_seq):
    """
    Modifies a DNA sequence to disrupt G-quadruplex structures.
    
    This function searches for occurrences of the "GGG" triplet in the DNA sequence
    and replaces them with synonymous codons. The aim is to disrupt G-quadruplex
    structures without altering the amino acid sequence encoded by the DNA.
    
    Args:
        DNA_seq (str): The DNA sequence to be modified.

    Returns:
        str: The modified DNA sequence with G-quadruplex structures disrupted.
    """

    # Split the DNA sequence into codons (groups of 3 bases)
    codon_seq = [DNA_seq[i:i+3] for i in range(0, len(DNA_seq), 3)]

    # Initialize a counter to limit the number of modification attempts
    attempt_count = 0

    # Loop until no more "GGG" triplets are found or the attempt limit is reached
    while True:
        attempt_count += 1

        # Break the loop if the number of attempts exceeds 100
        if attempt_count > 100:
            break

        # Search for the "GGG" triplet in the sequence
        pos = DNA_seq.find("GGG")
        
        # If "GGG" is not found, break the loop
        if pos == -1:
            break
        else:
            # Check the position of "GGG" in relation to the codon boundaries
            if pos % 3 == 0:  # "GGG" forms a complete codon
                # Replace the "GGG" codon with a random synonymous codon
                codon_seq[pos//3] = random.choice(dt.codon_synonyms["GGG"])
            elif pos % 3 == 1:  # "GGG" spans two codons, starting with the second base
                # Get the preceding codon
                codon1 = codon_seq[(pos-1)//3]
                # Replace the preceding codon with a random synonymous codon, if available
                if len(dt.codon_synonyms[codon1]) > 0:
                    codon_seq[(pos-1)//3] = random.choice(dt.codon_synonyms[codon1])
            elif pos % 3 == 2:  # "GGG" spans two codons, starting with the third base
                # Get the preceding codon
                codon1 = codon_seq[(pos-2)//3]
                # Replace the preceding codon with a random synonymous codon, if available
                if len(dt.codon_synonyms[codon1]) > 0:
                    codon_seq[(pos-2)//3] = random.choice(dt.codon_synonyms[codon1])

    # Join the modified codons into a single string and return the modified sequence
    return ''.join(codon_seq)


################################
# Remove 'AAAAA', 'CCCCC', etc #
################################


def strip_mononucleotide_tracts(DNA_seq):
    """
    Modifies a DNA sequence to disrupt mononucleotide tracts.
    
    This function identifies and alters tracts of five or more identical nucleotides 
    in a DNA sequence. It replaces codons adjacent to or containing these tracts 
    with synonymous codons, aiming to disrupt the tracts while preserving the 
    encoded amino acid sequence.
    
    Args:
        DNA_seq (str): The DNA sequence to be modified.

    Returns:
        str: The modified DNA sequence with mononucleotide tracts disrupted.
    """

    # Split the DNA sequence into codons (groups of 3 bases)
    codon_seq = [DNA_seq[i:i+3] for i in range(0, len(DNA_seq), 3)]

    # Initialize a counter to limit the number of modification attempts
    attempt_count = 0

    # Loop to identify and modify mononucleotide tracts
    while True:
        attempt_count += 1

        # Break the loop if the number of attempts exceeds 1000
        if attempt_count > 1000:
            break

        # Search for mononucleotide tracts of each base type
        pos = -1  # Initialize position variable
        for base in ["A", "C", "G", "T"]:
            tract = base * 5  # Create a string of 5 identical bases
            pos = DNA_seq.find(tract)
            if pos != -1:
                break

        # Break the loop if no mononucleotide tract is found
        if pos == -1:
            break

        # Identify and replace the codons adjacent to or containing the tract
        if pos % 3 == 0:  # Tract starts at the beginning of a codon
            codon1 = codon_seq[pos//3]
            if len(dt.codon_synonyms[codon1]) > 0:
                codon_seq[pos//3] = random.choice(dt.codon_synonyms[codon1])
        elif pos % 3 == 1:  # Tract spans two codons, starting with the second base
            for offset in [-1, 2]:
                codon_index = (pos + offset) // 3
                codon = codon_seq[codon_index]
                if len(dt.codon_synonyms[codon]) > 0:
                    codon_seq[codon_index] = random.choice(dt.codon_synonyms[codon])
        elif pos % 3 == 2:  # Tract spans two codons, starting with the third base
            for offset in [-2, 1]:
                codon_index = (pos + offset) // 3
                codon = codon_seq[codon_index]
                if len(dt.codon_synonyms[codon]) > 0:
                    codon_seq[codon_index] = random.choice(dt.codon_synonyms[codon])

    # Join the modified codons into a single string and return
    return ''.join(codon_seq)


##############################
# find+alter repeated 8-mers #
##############################

def strip_microhomology(DNA_seq):
    """
    Modifies a DNA sequence to disrupt microhomology.
    
    This function identifies and alters repeated 8-base sequences (8-mers)
    within the DNA sequence. When an 8-mer occurs more than once, the function
    replaces adjacent codons with synonymous codons to disrupt the microhomology
    while preserving the encoded amino acid sequence.
    
    Args:
        DNA_seq (str): The DNA sequence to be modified.

    Returns:
        str: The modified DNA sequence with microhomology disrupted.
    """

    # Split the DNA sequence into codons (groups of 3 bases)
    codon_seq = [DNA_seq[i:i+3] for i in range(0, len(DNA_seq), 3)]

    # Create a library of all 8-mers in the sequence
    kmer_lib = {}
    for i in range(len(DNA_seq)-8):
        kmer = DNA_seq[i:i+8]
        kmer_lib[kmer] = kmer_lib.get(kmer, 0) + 1

    # Find and modify repeated 8-mers
    for kmer in kmer_lib:
        # Skip 8-mers that occur only once
        if kmer_lib[kmer] == 1:
            continue

        # Find the position of the first occurrence of the repeated 8-mer
        pos = DNA_seq.find(kmer)

        # Identify and replace codons adjacent to or containing the 8-mer
        if pos % 3 == 0:  # 8-mer starts at the beginning of a codon
            for offset in [0, 1]:
                codon_index = (pos // 3) + offset
                codon = codon_seq[codon_index]
                if len(dt.codon_synonyms[codon]) > 0:
                    codon_seq[codon_index] = random.choice(dt.codon_synonyms[codon])

        elif pos % 3 == 1:  # 8-mer spans two codons, starting with the second base
            for offset in [-1, 2, 3]:
                codon_index = (pos + offset) // 3
                codon = codon_seq[codon_index]
                if len(dt.codon_synonyms[codon]) > 0:
                    codon_seq[codon_index] = random.choice(dt.codon_synonyms[codon])

        elif pos % 3 == 2:  # 8-mer spans two codons, starting with the third base
            for offset in [-2, 1, 2]:
                codon_index = (pos + offset) // 3
                codon = codon_seq[codon_index]
                if len(dt.codon_synonyms[codon]) > 0:
                    codon_seq[codon_index] = random.choice(dt.codon_synonyms[codon])

        # Decrement the count for the modified 8-mer
        kmer_lib[kmer] -= 1

    # Join the modified codons into a single string and return
    return ''.join(codon_seq)

#######################
# Improve G-C content #
#######################


def fix_GC(
        DNA_seq,
        lim_high: float=0.6,
        lim_low: float=0.45
    ):
    """
    Adjusts the GC content of a DNA sequence within specified limits.

    This function modifies the DNA sequence by replacing codons with synonymous 
    codons to ensure the GC content is within the specified range. It aims to 
    decrease GC content if it's above lim_high, and increase it if below lim_low, 
    while preserving the amino acid sequence encoded by the DNA.

    Args:
        DNA_seq (str): The DNA sequence to be modified.
        lim_high (float): Upper limit of desired GC content.
        lim_low (float): Lower limit of desired GC content.

    Returns:
        str: Modified DNA sequence with adjusted GC content.
    """
    
    # Calculate initial GC content
    GC_init = (DNA_seq.count("G") + DNA_seq.count("C")) / len(DNA_seq)

    # Check if GC content adjustment is needed
    if lim_low <= GC_init <= lim_high:
        return DNA_seq

    # Split sequence into codons
    codon_seq = [DNA_seq[i:i+3] for i in range(0, len(DNA_seq), 3)]

    # Adjust GC content if it's too high or too low
    if GC_init > lim_high:
        target = 'decrease'
    elif GC_init < lim_low:
        target = 'increase'

    attempt_count = 0
    while True:
        attempt_count += 1

        # Terminate if it's impossible to adjust GC content
        if attempt_count > 1000:
            sys.stderr.write("Failed to fix GC to within specified limits\n")
            sys.exit()

        # Select a random codon for potential replacement
        pos = random.randrange(0, len(codon_seq))
        codon_init = codon_seq[pos]

        # Shuffle the list of synonymous codons
        synonymous_codons = dt.codon_synonyms[codon_init][:]
        random.shuffle(synonymous_codons)

        # Attempt to find a suitable synonymous codon
        for codon_new in synonymous_codons:
            GC_new = (codon_new.count("G") + codon_new.count("C")) / 3
            GC_init_codon = (codon_init.count("G") + codon_init.count("C")) / 3

            if (target == 'decrease' and GC_new < GC_init_codon) or \
               (target == 'increase' and GC_new > GC_init_codon):
                codon_seq[pos] = codon_new
                break

        # Reconstruct the DNA sequence and check if GC content is within limits
        modified_seq = ''.join(codon_seq)
        GC_modified = (modified_seq.count("G") + modified_seq.count("C")) / len(modified_seq)
        
        if (target == 'decrease' and GC_modified < lim_high) or \
           (target == 'increase' and GC_modified > lim_low):
            return modified_seq
