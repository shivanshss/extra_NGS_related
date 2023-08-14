#Given: An array A for which A[k] represents the proportion of homozygous recessive individuals for the k-th Mendelian factor in a diploid population. Assume that the population is in genetic equilibrium for all factors.

#Return: An array B having the same length as A in which B[k] represents the probability that a randomly selected individual carries at least one copy of the recessive allele for the k-th factor.

def calculate_probabilities(A):
    B = []
    for p in A:
        q = 1 - p  # Probability of being homozygous dominant
        B.append(1 - q**2)  # Probability of carrying at least one recessive allele
    return B

# Sample input
A = [0.1, 0.25, 0.5]

# Calculate probabilities
B = calculate_probabilities(A)

# Print the result
print(" ".join(map(str, B)))

