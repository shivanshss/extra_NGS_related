#Given: Two protein strings s and t in FASTA format (each of length at most 100 aa).

#Return: The maximum alignment score between s and t, followed by two augmented strings s′ and t′ representing an optimal alignment of s and t

#Use:
#The BLOSUM62 scoring matrix.
#Gap opening penalty equal to 11.
#Gap extension penalty equal to 1.


from Bio import SeqIO
from Bio.Align import substitution_matrices
from Bio import pairwise2

def max_alignment_score(s, t, gap_open, gap_extend):
    # Define the BLOSUM62 scoring matrix
    matrix = substitution_matrices.load("BLOSUM62")

    # Perform global sequence alignment with affine gap penalty
    alignments = pairwise2.align.globaldx(s, t, matrix, gap_open, gap_extend, penalize_extend_when_opening=True)

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

    # Define the gap penalties
    gap_open = 11
    gap_extend = 1

    # Calculate the maximum alignment score
    max_score = max_alignment_score(s, t, gap_open, gap_extend)
    print(max_score)

if __name__ == "__main__":
    main()

