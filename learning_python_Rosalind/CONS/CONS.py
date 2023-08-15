from Bio import Align
import numpy as np

def align_sequences(dna_strings):
    aligner = Align.PairwiseAligner()
    alignments = aligner.align(dna_strings)
    aligned_sequences = [alignment[0] for alignment in alignments]
    return aligned_sequences

def build_profile_matrix(aligned_sequences):
    profile = {'A': [], 'C': [], 'G': [], 'T': []}
    alignment_length = len(aligned_sequences[0])
    
    for nucleotide in 'ACGT':
        for i in range(alignment_length):
            count = sum(seq[i] == nucleotide for seq in aligned_sequences)
            profile[nucleotide].append(count)
    
    return profile

def generate_consensus_string(profile):
    consensus = ""
    alignment_length = len(profile['A'])
    
    for i in range(alignment_length):
        max_count = max(profile[nucleotide][i] for nucleotide in 'ACGT')
        consensus_nucleotide = [nucleotide for nucleotide in 'ACGT' if profile[nucleotide][i] == max_count][0]
        consensus += consensus_nucleotide
    
    return consensus

def main(input_file_path):
    with open(input_file_path, 'r') as input_file:
        dna_strings = []
        current_string = ""
        
        for line in input_file:
            line = line.strip()
            if line.startswith('>'):
                if current_string:
                    dna_strings.append(current_string)
                current_string = ""
            else:
                current_string += line
        dna_strings.append(current_string)
        
    aligned_sequences = align_sequences(dna_strings)
    profile = build_profile_matrix(aligned_sequences)
    consensus = generate_consensus_string(profile)
    
    return consensus, profile

# Example usage
input_file_path = "input.fasta"  # Replace with your input file path
consensus, profile = main(input_file_path)

# Print the results
print("Consensus:")
print(consensus)
print("\nProfile:")
for nucleotide in 'ACGT':
    counts = " ".join(map(str, profile[nucleotide]))
    print(f"{nucleotide}: {counts}")

