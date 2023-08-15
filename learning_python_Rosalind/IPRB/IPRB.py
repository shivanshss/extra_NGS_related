#Given: Three positive integers k, m, and n, representing a population containing k+m+n organisms: k individuals are homozygous dominant for a factor, m are heterozygous, and n are homozygous recessive.

#Return: The probability that two randomly selected mating organisms will produce an individual possessing a dominant allele (and thus displaying the dominant phenotype). Assume that any two organisms can mate.


import pandas as pd

def dominant_probability(k, m, n):
    total_population = k + m + n

    # Probabilities of different genotype combinations
    prob_kk = (k / total_population) * ((k - 1) / (total_population - 1))
    prob_mm = (m / total_population) * ((m - 1) / (total_population - 1)) * 0.75
    prob_km = (k / total_population) * (m / (total_population - 1)) * 0.5
    prob_mk = (m / total_population) * (k / (total_population - 1)) * 0.5
    prob_kn = (k / total_population) * (n / (total_population - 1))
    prob_nk = (n / total_population) * (k / (total_population - 1))
    prob_mn = (m / total_population) * (n / (total_population - 1)) * 0.5
    prob_nm = (n / total_population) * (m / (total_population - 1)) * 0.5
    prob_nn = (n / total_population) * ((n - 1) / (total_population - 1)) * 0


    dominant_prob = prob_kk + prob_mm + prob_km + prob_mk + prob_kn + prob_nk + prob_mn + prob_nm + prob_nn
    return dominant_prob

# Read the input file into a pandas DataFrame
input_file_path = "input.txt"
df = pd.read_csv(input_file_path, sep='\s+', names=['k', 'm', 'n'])

# Calculate dominant allele probability for each set of values
df['probability'] = df.apply(lambda row: dominant_probability(row['k'], row['m'], row['n']), axis=1)

# Save the results to an output file
output_file_path = "output.txt"
df.to_csv(output_file_path, sep='\t', index=False)

print(df)

