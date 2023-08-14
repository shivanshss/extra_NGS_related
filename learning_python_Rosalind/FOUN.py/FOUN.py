#Given: Two positive integers N and m, followed by an array A containing k integers between 0 and 2N. A[j] represents the number of recessive alleles for the j-th factor in a population of N diploid individuals.

#Return: An m×k matrix B for which Bi,j represents the common logarithm of the probability that after i generations, no copies of the recessive allele for the j-th factor will remain in the population. Apply the Wright-Fisher model.

import math

def calculate_probability(N, m, A):
    B = []
    for i in range(m):
        row = []
        for j in range(len(A)):
            p = A[j] / (2 * N)  # Probability of selecting a recessive allele
            probability = (1 - p) ** i
            row.append(math.log10(probability))
        B.append(row)
    return B

# Sample input
N = 4
m = 3
A = [0, 1, 2]

# Calculate the matrix B
result_matrix = calculate_probability(N, m, A)

# Print the result
for row in result_matrix:
    print(" ".join(map(str, row)))

