def calculate_expected_offspring(genotypes):
    probabilities = [1, 1, 1, 0.75, 0.5, 0]  # Corresponding probabilities for each genotype pairing
    expected_offspring = sum(prob * count * 2 for prob, count in zip(probabilities, genotypes))
    return expected_offspring

def main(input_file_path):
    with open(input_file_path, 'r') as input_file:
        genotypes = list(map(int, input_file.readline().strip().split()))
    
    expected_offspring = calculate_expected_offspring(genotypes)
    return expected_offspring

# Example usage
input_file_path = "input.txt"  # Replace with your input file path
expected_offspring = main(input_file_path)
print(expected_offspring)

