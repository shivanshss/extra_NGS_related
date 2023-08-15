#Given: A DNA string s of length at most 1000 nt.

#Return: Four integers (separated by spaces) counting the respective number of times that the symbols 'A', 'C', 'G', and 'T' occur in s


def count_symbols(dna_string):
    symbol_count = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
    
    for symbol in dna_string:
        if symbol in symbol_count:
            symbol_count[symbol] += 1
            
    return symbol_count['A'], symbol_count['C'], symbol_count['G'], symbol_count['T']

def process_fasta_file(file_path):
    with open(file_path, 'r') as fasta_file:
        sequences = []
        current_sequence = ""
        
        for line in fasta_file:
            line = line.strip()
            if line.startswith('>'):
                if current_sequence:
                    sequences.append(current_sequence)
                current_sequence = ""
            else:
                current_sequence += line
        
        if current_sequence:
            sequences.append(current_sequence)
    
    for idx, sequence in enumerate(sequences, start=1):
        counts = count_symbols(sequence)
        print(f"Sequence {idx}: {counts[0]} {counts[1]} {counts[2]} {counts[3]}")

# Example usage
fasta_file_path = "input.fasta"  # Replace with the path to your FASTA file
process_fasta_file(fasta_file_path)

