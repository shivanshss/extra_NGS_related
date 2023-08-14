def transition_transversion_ratio(s1, s2):
    transitions = 0
    transversions = 0

    for base1, base2 in zip(s1, s2):
        if base1 != base2:
            if (base1 in 'AG' and base2 in 'AG') or (base1 in 'CT' and base2 in 'CT'):
                transitions += 1
            else:
                transversions += 1

    return transitions / transversions

def main():
    # Read input sequences from FASTA format
    sequences = []
    with open("input.fasta", "r") as fasta_file:
        for line in fasta_file:
            if not line.startswith(">"):
                sequences.append(line.strip())

    s1 = sequences[0]  # First DNA sequence
    s2 = sequences[1]  # Second DNA sequence

    ratio = transition_transversion_ratio(s1, s2)
    print(ratio)

if __name__ == "__main__":
    main()

