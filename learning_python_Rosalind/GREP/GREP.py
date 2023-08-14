#Given: A list Sk+1 of error-free DNA (k+1)-mers (k≤5) taken from the same strand of a circular chromosome (of length ≤50).

#Return: All circular strings assembled by complete cycles in the de Bruijn graph Bk of Sk+1. The strings may be given in any order, but each one should begin with the first (k+1)-mer provided in the input.

def construct_circular_strings(kmers):
    k = len(kmers[0]) - 1
    kmer_dict = {}
    for kmer in kmers:
        prefix = kmer[:-1]
        suffix = kmer[1:]
        if prefix in kmer_dict:
            kmer_dict[prefix].append(suffix)
        else:
            kmer_dict[prefix] = [suffix]

    circular_strings = []
    for start_kmer in kmer_dict.keys():
        current_kmer = start_kmer
        circular_string = current_kmer
        while True:
            next_suffixes = kmer_dict[current_kmer]
            if len(next_suffixes) != 1:
                break
            next_suffix = next_suffixes[0]
            circular_string += next_suffix[-1]
            current_kmer = next_suffix
            if current_kmer == start_kmer:
                circular_strings.append(circular_string)
                break

    return circular_strings

def main():
    kmers = []
    while True:
        try:
            line = input().strip()
            if line:
                kmers.append(line)
        except EOFError:
            break

    # Construct circular strings
    circular_strings = construct_circular_strings(kmers)
    for circular_string in circular_strings:
        print(circular_string)

if __name__ == "__main__":
    main()

