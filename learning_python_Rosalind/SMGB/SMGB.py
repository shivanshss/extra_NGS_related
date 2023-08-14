#Given: Two DNA strings s and t in FASTA format, each having length at most 10 kbp.

#Return: The maximum semiglobal alignment score of s and t, followed by an alignment of s and t achieving this maximum score. Use an alignment score in which matching symbols count +1, substitutions count -1, and there is a linear gap penalty of 1. If multiple optimal alignments exist, then you may return any one.

from Bio import SeqIO
from Bio import pairwise2

def semiglobal_alignment(s, t):
    # Define the scoring scheme
    match_score = 1
    mismatch_score = -1
    gap_penalty = -1

    # Perform semiglobal alignment
    alignments = pairwise2.align.globalds(s, t, match_score, mismatch_score, gap_penalty, gap_penalty, penalize_extend_when_opening=True)

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

    # Calculate the maximum semiglobal alignment score
    score, aligned_s, aligned_t = semiglobal_alignment(s, t)
    print(score)
    print(aligned_s)
    print(aligned_t)

if __name__ == "__main__":
    main()

