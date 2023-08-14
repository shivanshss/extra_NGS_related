#Given: Two DNA strings s and t in FASTA format that share some short inexact repeat r of 32-40 bp. By "inexact" we mean that r may appear with slight modifications (each repeat differ by ≤ 3 changes/indels).

# Return: The total number of occurrences of r as a substring of s, followed by the total number of occurrences of r as a substring of t.


from Bio import SeqIO
from Bio.SubsMat import MatrixInfo as matlist

def find_occurrences(seq, repeat, max_changes=3):
    occurrences = 0
    repeat_length = len(repeat)
    for i in range(len(seq) - repeat_length + 1):
        subseq = seq[i:i+repeat_length]
        num_changes = sum(a != b for a, b in zip(subseq, repeat))
        if num_changes <= max_changes:
            occurrences += 1
    return occurrences

def main():
    # Read input sequences from FASTA format
    sequences = []
    with open("input.fasta", "r") as fasta_file:
        for record in SeqIO.parse(fasta_file, "fasta"):
            sequences.append(str(record.seq))

    # Extract sequences
    s = sequences[0]
    t = sequences[1]

    # Define the inexact repeat and its length range
    repeat_min_length = 32
    repeat_max_length = 40

    # Find the occurrences of the repeat in both sequences
    repeat_occurrences_s = 0
    repeat_occurrences_t = 0
    for repeat_length in range(repeat_min_length, repeat_max_length + 1):
        repeat = s[:repeat_length]  # Extract the repeat from sequence s
        repeat_occurrences_s += find_occurrences(s, repeat)
        repeat_occurrences_t += find_occurrences(t, repeat)

    # Print the results
    print(repeat_occurrences_s, repeat_occurrences_t)

if __name__ == "__main__":
    main()

