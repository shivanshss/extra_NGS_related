#Given: Two protein strings s and t in FASTA format, each of length at most 1000 aa.

#Return: The total number of optimal alignments of s and t with respect to edit alignment score, modulo 134,217,727 (227-1).


def count_optimal_alignments(s, t):
    m = len(s)
    n = len(t)

    # Create a 2D table to store alignment counts
    counts = [[0] * (n + 1) for _ in range(m + 1)]

    # Initialize the first row and column
    for i in range(m + 1):
        counts[i][0] = 1
    for j in range(n + 1):
        counts[0][j] = 1

    # Fill in the table using dynamic programming
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            if s[i - 1] == t[j - 1]:
                counts[i][j] = counts[i - 1][j - 1]
            else:
                counts[i][j] = (counts[i - 1][j] + counts[i][j - 1]) % 134217727

    return counts[m][n]

def main():
    # Read input sequences from FASTA format
    sequences = []
    with open("input.fasta", "r") as fasta_file:
        lines = fasta_file.readlines()
        sequences = [line.strip() for line in lines[1:]]

    s = sequences[0]  # First protein string
    t = sequences[1]  # Second protein string

    # Count the total number of optimal alignments
    num_alignments = count_optimal_alignments(s, t)
    print(num_alignments)

if __name__ == "__main__":
    main()

