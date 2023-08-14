#Given: A collection of up to 1000 (possibly repeating) DNA strings of equal length (not exceeding 50 bp) corresponding to a set S of (k+1) -mers.

#Return: The adjacency list corresponding to the de Bruijn graph corresponding to S∪Src.

def construct_de_bruijn_graph(strings, k):
    de_bruijn_graph = {}

    for s in strings:
        for i in range(len(s) - k + 1):
            kmer = s[i:i + k]
            if kmer not in de_bruijn_graph:
                de_bruijn_graph[kmer] = []
            if i + k < len(s):
                next_kmer = s[i + 1:i + k + 1]
                de_bruijn_graph[kmer].append(next_kmer)

    return de_bruijn_graph

def main():
    k = int(input())  # Length of k-mers
    strings = []
    while True:
        try:
            line = input().strip()
            if line:
                strings.append(line)
        except EOFError:
            break

    # Construct the de Bruijn graph
    de_bruijn_graph = construct_de_bruijn_graph(strings, k)

    # Print the adjacency list
    for node, neighbors in de_bruijn_graph.items():
        for neighbor in neighbors:
            print(f"({node}, {neighbor})")

if __name__ == "__main__":
    main()

