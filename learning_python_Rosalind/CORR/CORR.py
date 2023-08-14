#Given: A collection of up to 1000 reads of equal length (at most 50 bp) in FASTA format. Some of these reads were generated with a single-nucleotide error. For each read s in the dataset, one of the following applies:
#s was correctly sequenced and appears in the dataset at least twice (possibly as a reverse complement);
#s is incorrect, it appears in the dataset exactly once, and its Hamming distance is 1 with respect to exactly one correct read in the dataset (or its reverse complement).

#Return: A list of all corrections in the form "[old read]->[new read]". (Each correction must be a single symbol substitution, and you may return the corrections in any order.)

def find_sequencing_errors(reads):
    error_corrections = []

    # Create a dictionary to store read frequencies
    read_counts = {}
    for read in reads:
        if read in read_counts:
            read_counts[read] += 1
        else:
            read_counts[read] = 1

    # Find sequencing errors and their corrections
    for read in reads:
        if read_counts[read] == 1:
            for correct_read in read_counts:
                if correct_read != read and hamming_distance(correct_read, read) == 1:
                    error_corrections.append(f"{read}->{correct_read}")
                    break

    return error_corrections
    
def hamming_distance(s1, s2):
    if len(s1) > len(s2):
        s1, s2 = s2, s1  # Swap s1 and s2 if s1 is longer than s2

    dH_s1s2 = sum(1 for x, y in zip(s1, s2) if x != y)  # Calculate dH(s1, s2)

    dH_s2s1 = dH_s1s2 + len(s2) - len(s1)  # Calculate dH(s2, s1)
    
    return dH_s1s2, dH_s2s1
    
def hamming_distance(s1, s2):
    return sum(c1 != c2 for c1, c2 in zip(s1, s2))

def main():
    # Read input sequences from FASTA format
    sequences = []
    with open("input.fasta", "r") as fasta_file:
        lines = fasta_file.readlines()
        sequences = [line.strip() for line in lines[1:]]

    # Calculate the error corrections
    corrections = find_sequencing_errors(sequences)
    for correction in corrections:
        print(correction)

if __name__ == "__main__":
    main()

