def rabbit_pairs_after_n_months(n, m):
    dp = [0] * (m - 1) + [1]  # Initialize with 1 pair in the last generation
    for i in range(1, n):
        new_pairs = sum(dp[:-1])  # Sum of all pairs except the oldest generation
        dp[:-1] = dp[1:]  # Shift generations
        dp[-1] = new_pairs  # Add the new generation
    
    total_pairs = sum(dp)
    return total_pairs

def process_dataframe_file(input_file_path, output_file_path):
    import pandas as pd
    
    df = pd.read_csv(input_file_path)
    results = []

    for index, row in df.iterrows():
        n = row['n']
        m = row['m']
        total_pairs = rabbit_pairs_after_n_months(n, m)
        results.append((n, m, total_pairs))
    
    with open(output_file_path, 'w') as output_file:
        for result in results:
            output_file.write(f"{result[0]} {result[1]} {result[2]}\n")

# Example usage
input_file_path = "input.csv"   # Replace with your input CSV file path
output_file_path = "output.txt"  # Replace with your desired output file path
process_dataframe_file(input_file_path, output_file_path)

