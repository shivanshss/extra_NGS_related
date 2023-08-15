#Given: A DNA string t having length at most 1000 nt.

#Return: The transcribed RNA string of t

def transcribe_dna_to_rna(dna_string):
    rna_string = dna_string.replace('T', 'U')
    return rna_string

def process_multi_fasta_file(file_path):
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
        rna_sequence = transcribe_dna_to_rna(sequence)
        print(f"RNA Sequence {idx}: {rna_sequence}")

# Example usage
fasta_file_path = "input.fasta"
process_multi_fasta_file(fasta_file_path)

