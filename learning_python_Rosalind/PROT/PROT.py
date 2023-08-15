#Given: An RNA string s corresponding to a strand of mRNA (of length at most 10 kbp).

#Return: The protein string encoded by s


def rna_to_protein(rna_string):
    codon_table = {
    "UUU": "F", "UUC": "F", "UUA": "L", "UUG": "L",
    "CUU": "L", "CUC": "L", "CUA": "L", "CUG": "L",
    "AUU": "I", "AUC": "I", "AUA": "I", "AUG": "M",
    "GUU": "V", "GUC": "V", "GUA": "V", "GUG": "V",
    "UCU": "S", "UCC": "S", "UCA": "S", "UCG": "S",
    "CCU": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACU": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCU": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "UAU": "Y", "UAC": "Y", "UAA": "*", "UAG": "*",
    "CAU": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAU": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAU": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "UGU": "C", "UGC": "C", "UGA": "*", "UGG": "W",
    "CGU": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGU": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGU": "G", "GGC": "G", "GGA": "G", "GGG": "G"
    }


    protein = ""
    for i in range(0, len(rna_string), 3):
        codon = rna_string[i:i+3]
        amino_acid = codon_table.get(codon, "")
        if amino_acid == "*":
            break  # Stop at the stop codon
        protein += amino_acid

    return protein

def process_multi_fasta_file(file_path):
    with open(file_path, 'r') as fasta_file:
        sequences = {}
        current_label = ""
        
        for line in fasta_file:
            line = line.strip()
            if line.startswith('>'):
                current_label = line[1:]
                sequences[current_label] = ""
            else:
                sequences[current_label] += line
    
    for label, rna_string in sequences.items():
        protein_string = rna_to_protein(rna_string)
        print(f">{label}")
        print(protein_string)

# Example usage
fasta_file_path = "input.fasta"  # Replace with the path to your multi-FASTA file
process_multi_fasta_file(fasta_file_path)

