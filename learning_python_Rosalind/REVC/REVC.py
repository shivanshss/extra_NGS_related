#Given: A DNA string s of length at most 1000 bp.

#Return: The reverse complement sc of s.


def reverse_complement(dna_string):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    reverse_comp_string = ''.join(complement[symbol] for symbol in reversed(dna_string))
    return reverse_comp_string

def process_multi_fasta_file(input_file_path, output_file_path):
    with open(input_file_path, 'r') as input_fasta, open(output_file_path, 'w') as output_fasta:
        sequences = []
        current_sequence = ""
        
        for line in input_fasta:
            line = line.strip()
            if line.startswith('>'):
                if current_sequence:
                    sequences.append(current_sequence)
                current_sequence = ""
                output_fasta.write(line + "\n")
            else:
                current_sequence += line
        
        if current_sequence:
            sequences.append(current_sequence)
    
        for sequence in sequences:
            reverse_comp_sequence = reverse_complement(sequence)
            output_fasta.write(reverse_comp_sequence + "\n")

# Example usage
input_fasta_path = "input.fasta"
output_fasta_path = "output.fasta"
process_multi_fasta_file(input_fasta_path, output_fasta_path)

