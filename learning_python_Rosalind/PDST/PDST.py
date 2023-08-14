#Given: A collection of n (n≤10) DNA strings s1,…,sn of equal length (at most 1 kbp). Strings are given in FASTA format.

#Return: The matrix D corresponding to the p-distance dp on the given strings. As always, note that your answer is allowed an absolute error of 0.001.

from Bio import SeqIO

def p_distance(s1, s2):
    differences = sum(a != b for a, b in zip(s1, s2))
    return differences / len(s1)

def calculate_distance_matrix(sequences):
    n = len(sequences)
    distance_matrix = [[0.0] * n for _ in range(n)]

    for i in range(n):
        for j in range(i + 1, n):
            distance = p_distance(sequences[i], sequences[j])
            distance_matrix[i][j] = distance
            distance_matrix[j][i] = distance

    return distance_matrix

def main():
    # Read input sequences from FASTA format
    sequences = list(SeqIO.parse("input.fasta", "fasta"))

    # Calculate the p-distance matrix
    distance_matrix = calculate_distance_matrix(sequences)

    # Print the distance matrix with formatting
    for row in distance_matrix:
        print(" ".join([f"{dist:.5f}" for dist in row]))

if __name__ == "__main__":
    main()

