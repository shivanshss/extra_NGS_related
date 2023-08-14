#Given: Two DNA strings s1 and s2 of equal length (at most 1 kbp).

#Return: The transition/transversion ratio R(s1,s2)


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
    with open("input.fasta", "r") as fasta_file:
        sequences = fasta_file.read().splitlines()

    s1 = sequences[1]  # First DNA sequence
    s2 = sequences[3]  # Second DNA sequence

    ratio = transition_transversion_ratio(s1, s2)
    print(ratio)

if __name__ == "__main__":
    main()

