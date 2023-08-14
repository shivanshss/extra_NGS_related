#Given: Two DNA strings s and t in FASTA format, each having length at most 10 kbp.

#Return: The score of an optimal overlap alignment of s and t, followed by an alignment of a suffix s′ of s and a prefix t′ of t achieving this optimal score. Use an alignment score in which matching symbols count +1, substitutions count -2, and there is a linear gap penalty of 2. If multiple optimal alignments exist, then you may return any one.

from Bio import SeqIO
from Bio import pairwise2

def optimal_overlap_alignment(s, t):
    # Define the scoring scheme
    match_score = 1
    mismatch_score = -2
    gap_open_penalty = -2
    gap_extend_penalty = -2

    # Perform local sequence alignment
    alignments = pairwise2.align.localds(s, t, match_score, mismatch_score, gap_open_penalty, gap_extend_penalty)

    # Get the alignment with the maximum score
    max_alignment = max(alignments, key=lambda x: x.score)

    return max_alignment.score, max_alignment.seqA, max_alignment.seqB

def main():
    # Read input sequences from FASTA format
    sequences = []
    with open("input.fasta", "r") as fasta_file:
        lines = fasta_file.readlines()
        sequences = [line.strip() for line in lines[1:]]

    s = sequences[0]  # First DNA string
    t = sequences[1]  # Second DNA string

    # Calculate the optimal overlap alignment
    score, aligned_s, aligned_t = optimal_overlap_alignment(s, t)
    print(score)
    print(aligned_s)
    print(aligned_t)

if __name__ == "__main__":
    main()

