#install Biopython package 
pip install Biopython


#import pakcages 
from Bio import SeqIO
from Bio.Seq import Seq
from io import StringIO

# Read fasta file here i am using colab so path is look liike this but it need to change with respect to location of file 
input_file = '/content/influenza.fna'

# To extract 10 sequences from multiple sequence fasta file 
extracted_sequences, remaining_sequences = extract_sequences(input_file, num_sequences=10)

def extract_sequences(input_file, num_sequences=10):
    records = list(SeqIO.parse(input_file, "fasta"))
    extracted_sequences = records[:num_sequences]
    remaining_sequences = records[num_sequences:]
    return extracted_sequences, remaining_sequences

def extract_sequences(input_file, num_sequences=10):
    records = list(SeqIO.parse(input_file, "fasta"))
    extracted_sequences = records[:num_sequences]
    remaining_sequences = records[num_sequences:]
    return extracted_sequences, remaining_sequences


def reverse_complement(sequence):
    return str(Seq(sequence).reverse_complement())


modified_sequences = []
for seq_record in extracted_sequences:
    original_sequence = str(seq_record.seq)
    reversed_complement_sequence = reverse_complement(original_sequence)
    modified_sequence = reversed_complement_sequence[:-1] + 'A'  # Change the last nucleotide to 'A'
    modified_sequences.append(SeqIO.SeqRecord(Seq(modified_sequence), id=seq_record.id, description=seq_record.description))

#lets see output now 
output_prefix = 'influenza'

# Save the extracted sequences to a new file
output_file_path = f'{output_prefix}_extracted.fasta'
SeqIO.write(extracted_sequences, output_file_path, "fasta")

# Create BED file from extracted sequences
bed_content = create_bed_from_sequences(extracted_sequences)

# Print BED content
print(bed_content)

with open('influenza.bed', 'w') as bed_file:
    bed_file.write(bed_content)

# Insert modified sequences back into the remaining sequences
all_sequences = remaining_sequences + modified_sequences

# Print the modified sequences
for seq_record in all_sequences:
    print(f">{seq_record.id}\n{seq_record.seq}")



modified_output_file_path = f'{output_prefix}_modified.fasta'
SeqIO.write(all_sequences, modified_output_file_path, "fasta")



















