#Given: Positive integers N (N≤7), m (m≤2N), g (g≤6) and k (k≤2N).

#Return: The probability that in a population of N diploid individuals initially possessing m copies of a dominant allele, we will observe after g generations at least k copies of a recessive allele. Assume the Wright-Fisher model.

from math import comb

def calculate_probability(N, m, g, k):
    p = m / (2 * N)  # Probability of selecting a recessive allele
    probability = 0

    for r in range(k, g + 1):
        probability += comb(g, r) * (p ** r) * ((1 - p) ** (g - r))

    return probability

# Sample input
N = 4
m = 6
g = 2
k = 1

# Calculate the probability
result = calculate_probability(N, m, g, k)

# Print the result
print(result)

