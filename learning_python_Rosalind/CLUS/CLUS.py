#Given: Set of nucleotide strings in FASTA format.

#Return: ID of the string most different from the others.

from Bio import SeqIO
from Bio.Align import pairwise2

def calculate_average_identity(seq1, seq2):
    matches = 0
    total_positions = len(seq1)

    for base1, base2 in zip(seq1, seq2):
        if base1 == base2:
            matches += 1

    return matches / total_positions

def most_different_sequence(sequences):
    max_avg_identity = 0
    most_different_id = ""

    for i, seq1 in enumerate(sequences):
        avg_identity_sum = 0

        for j, seq2 in enumerate(sequences):
            if i != j:  # Skip comparing the same sequence
                avg_identity_sum += calculate_average_identity(seq1, seq2)

        avg_identity = avg_identity_sum / (len(sequences) - 1)

        if avg_identity > max_avg_identity:
            max_avg_identity = avg_identity
            most_different_id = seq1.id

    return most_different_id

def main():
    # Read input sequences from FASTA format
    sequences = list(SeqIO.parse("input.fasta", "fasta"))

    # Find the most different sequence
    most_different_id = most_different_sequence(sequences)
    print(most_different_id)

if __name__ == "__main__":
    main()

