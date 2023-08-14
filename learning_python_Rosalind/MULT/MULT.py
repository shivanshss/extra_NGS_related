#Given: A collection of four DNA strings of length at most 10 bp in FASTA format.

#Return: A multiple alignment of the strings having maximum score, where we score matched symbols 0 (including matched gap symbols) and all mismatched symbols -1 (thus incorporating a linear gap penalty of 1).


from Bio import SeqIO
from Bio.Align import MultipleSeqAlignment
from Bio import pairwise2

def score_function(x, y):
    if x == y:
        return 0
    else:
        return -1

def max_alignment_score(strings):
    alignments = pairwise2.align.multiple(strings, gap_char="-", score_only=True, score_function=score_function)
    return sum(alignments)

def main():
    # Read input sequences from FASTA format
    sequences = []
    with open("input.fasta", "r") as fasta_file:
        lines = fasta_file.readlines()
        sequences = [line.strip() for line in lines[1:]]

    # Calculate the multiple alignment score
    max_score = max_alignment_score(sequences)
    print(max_score)

if __name__ == "__main__":
    main()

