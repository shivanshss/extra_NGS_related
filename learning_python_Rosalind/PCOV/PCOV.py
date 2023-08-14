#Given: A collection of (error-free) DNA k-mers (k≤50) taken from the same strand of a circular chromosome. In this dataset, all k-mers from this strand of the chromosome are present, and their de Bruijn graph consists of exactly one simple cycle.

#Return: A cyclic superstring of minimal length containing the reads (thus corresponding to a candidate cyclic chromosome).

def construct_cyclic_superstring(kmers):
    # Construct a dictionary to store k-1 mers as keys and their suffix as values
    suffix_dict = {}
    for kmer in kmers:
        prefix = kmer[:-1]
        suffix = kmer[1:]
        if prefix in suffix_dict:
            suffix_dict[prefix].append(suffix)
        else:
            suffix_dict[prefix] = [suffix]

    # Construct the cyclic superstring
    current_node = next(iter(suffix_dict.keys()))
    cyclic_superstring = current_node
    while True:
        next_suffixes = suffix_dict[current_node]
        next_suffix = next_suffixes.pop()
        if not next_suffixes:
            del suffix_dict[current_node]
        cyclic_superstring += next_suffix[-1]
        if next_suffix in suffix_dict:
            current_node = next_suffix
        else:
            break

    return cyclic_superstring

def main():
    kmers = []
    while True:
        try:
            line = input().strip()
            if line:
                kmers.append(line)
        except EOFError:
            break

    # Construct the cyclic superstring
    cyclic_superstring = construct_cyclic_superstring(kmers)
    print(cyclic_superstring)

if __name__ == "__main__":
    main()

