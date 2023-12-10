# Minimalist Codon Optimization Translation and DNA Design Library Preparation for Proteins

## Purpose

This repository is dedicated to converting libraries of designed protein sequences, including mutant variants, into corresponding DNA sequences. A key feature is the use of species-optimized codons, facilitating effective translation and transcription processes.

## How to Run the Library

### Step 1: Reverse Translation and Transcription
If you're starting with a library of amino acid-based protein sequences, the first step involves reverse translation and transcription to convert these sequences into unoptimized DNA.

To execute this process, run the following shell script:

```
cd scripts
sh run_reverse_aa_translation.sh
```

*Note*: Variable descriptions (e.g., `input` for the input FASTA file) are provided within the `.sh` script.

### Step 2: DNA Sequence Optimization
After converting proteins into DNA sequences, the next step involves optimization. This process includes codon optimization for a specific species, GC content optimization, microhomology region removal, mononucleotide tract removal, and G-quadruplex stripping.

Run the following shell script to initiate this process:

```
cd scripts
sh run_main.sh
```

*Note*: Details about input variables can be found within the shell script.

## Acknowledgements

Hat tip to Bryan Phillips for his generous contribution in sharing his repository, which includes minimalist code for reverse amino acid translation and species-specific codon optimization.
