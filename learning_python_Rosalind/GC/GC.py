#Given: At most 10 DNA strings in FASTA format (of length at most 1 kbp each).

#Return: The ID of the string having the highest GC-content, followed by the GC-content of that string. Rosalind allows for a default error of 0.001 in all decimal answers unless otherwise stated; please see the note on absolute error below.




def calculate_gc_content(dna_string):
    gc_count = dna_string.count('G') + dna_string.count('C')
    total_bases = len(dna_string)
    gc_content = (gc_count / total_bases) * 100
    return gc_content

def process_multi_fasta_file(input_file_path, output_file_path):
    fasta_data = {}
    current_label = ""
    
    with open(input_file_path, 'r') as fasta_file:
        for line in fasta_file:
            line = line.strip()
            if line.startswith('>'):
                current_label = line[1:]
                fasta_data[current_label] = ""
            else:
                fasta_data[current_label] += line
    
    highest_gc_id, highest_gc_content = "", 0.0

    for label, dna_string in fasta_data.items():
        gc_content = calculate_gc_content(dna_string)
        if gc_content > highest_gc_content:
            highest_gc_content = gc_content
            highest_gc_id = label

    with open(output_file_path, 'w') as output_file:
        output_file.write(highest_gc_id + '\n')
        output_file.write("{:.6f}".format(highest_gc_content) + '\n')

# Example usage
input_fasta_path = "input.fasta"
output_file_path = "output.txt"
process_multi_fasta_file(input_fasta_path, output_file_path)

