# Handling-FASTA-
generate bed from selective fasta 


This script performs various operations on a FASTA file containing Influenza  DNA sequences.
It uses the Biopython library for sequence manipulation.

Functions:
1. extract_sequences(input_file, num_sequences=10):
    Extracts a specified number of sequences from the input FASTA file.
    Returns the extracted sequences and the remaining sequences.

2. create_bed_from_sequences(sequences):
    Generates BED format content from a list of sequence records.
    BED format represents genomic features with chromosomal coordinates.

3. reverse_complement(sequence):
    Computes the reverse complement of a given DNA sequence.

Main Script:
1. Define input and output file paths:
   'Influenza.fna' file ware taken as input and influenza_extracted.fasta, infulenza.bed and influenza_modified.fasta as output 

2. Extract 10 sequences from the input FASTA file:

3. Save the extracted sequences to a new FASTA file:

4. Create BED file from extracted sequences:

5. Print BED content and save to a file:


6. Reverse complement and make a single nucleotide change:
    - For each extracted sequence, the last nucleotide is changed to 'A'.

7. Insert modified sequences back into the remaining sequences:
 
8. Print modified sequences and save to a new FASTA file:
