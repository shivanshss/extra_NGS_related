#Given: Two protein strings s and t in FASTA format (each of length at most 1000 aa).

#Return: The maximum alignment score between s and t

#Use:
#The BLOSUM62 scoring matrix.
#Constant gap penalty equal to 5.


from Bio import SeqIO
from Bio.Align import substitution_matrices
from Bio import pairwise2

def max_alignment_score(s, t):
    # Define the BLOSUM62 scoring matrix
    matrix = substitution_matrices.load("BLOSUM62")

    # Define the constant gap penalty
    gap_penalty = -5

    # Perform global sequence alignment
    alignments = pairwise2.align.globalds(s, t, matrix, gap_penalty, gap_penalty)

    # Get the maximum alignment score from the alignments
    max_score = max(alignment.score for alignment in alignments)

    return max_score

def main():
    # Read input sequences from FASTA format
    sequences = []
    with open("input.fasta", "r") as fasta_file:
        lines = fasta_file.readlines()
        sequences = [line.strip() for line in lines[1:]]

    s = sequences[0]  # First protein string
    t = sequences[1]  # Second protein string

    # Calculate the maximum alignment score
    max_score = max_alignment_score(s, t)
    print(max_score)

if __name__ == "__main__":
    main()

