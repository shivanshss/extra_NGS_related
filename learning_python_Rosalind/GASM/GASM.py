#Given: A collection S of (error-free) reads of equal length (not exceeding 50 bp). In this dataset, for some positive integer k, the de Bruijn graph Bk on Sk+1∪Srck+1 consists of exactly two directed cycles.

#Return: A cyclic superstring of minimal length containing every read or its reverse complement.

def construct_cyclic_superstring(reads):
    # Construct a set of all k-mers and their reverse complements
    k = len(reads[0])
    kmers = set()
    for read in reads:
        kmers.add(read)
        kmers.add(reverse_complement(read))

    # Create a dictionary to store k-1 mers as keys and their suffix as values
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

def reverse_complement(seq):
    complement_dict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join(complement_dict[base] for base in reversed(seq))

def main():
    reads = []
    while True:
        try:
            line = input().strip()
            if line:
                reads.append(line)
        except EOFError:
            break

    # Construct the cyclic superstring
    cyclic_superstring = construct_cyclic_superstring(reads)
    print(cyclic_superstring)

if __name__ == "__main__":
    main()

