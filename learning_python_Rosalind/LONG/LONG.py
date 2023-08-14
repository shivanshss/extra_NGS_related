#Given: At most 50 DNA strings of approximately equal length, not exceeding 1 kbp, in FASTA format (which represent reads deriving from the same strand of a single linear chromosome). The dataset is guaranteed to satisfy the following condition: there exists a unique way to reconstruct the entire chromosome from these reads by gluing together pairs of reads that overlap by more than half their length.

#Return: A shortest superstring containing all the given strings (thus corresponding to a reconstructed chromosome).


def find_shortest_superstring(strings):
    while len(strings) > 1:
        max_overlap = 0
        best_i = 0
        best_j = 0
        best_merged = ""

        for i in range(len(strings)):
            for j in range(len(strings)):
                if i != j:
                    overlap = 0
                    for k in range(1, min(len(strings[i]), len(strings[j])) + 1):
                        if strings[i][-k:] == strings[j][:k]:
                            overlap = k

                    if overlap > max_overlap:
                        max_overlap = overlap
                        best_i = i
                        best_j = j
                        best_merged = strings[i] + strings[j][max_overlap:]

        strings.pop(best_i)
        strings.pop(best_j - 1)
        strings.append(best_merged)

    return strings[0]

def main():
    # Read input sequences from FASTA format
    sequences = []
    with open("input.fasta", "r") as fasta_file:
        lines = fasta_file.readlines()
        sequences = [line.strip() for line in lines[1:]]

    # Calculate the shortest superstring
    superstring = find_shortest_superstring(sequences)
    print(superstring)

if __name__ == "__main__":
    main()

