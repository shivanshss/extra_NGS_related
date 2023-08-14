#Given: A collection of at most 1000 DNA strings (whose combined length does not exceed 50 kbp).

#Return: N50 and N75 for this collection of strings.

def calculate_n_statistic(contigs, percent):
    total_length = sum(len(contig) for contig in contigs)
    threshold_length = total_length * percent // 100

    contigs.sort(key=len, reverse=True)
    cumulative_length = 0
    for contig in contigs:
        cumulative_length += len(contig)
        if cumulative_length >= threshold_length:
            return len(contig)

def main():
    contigs = []
    while True:
        try:
            line = input().strip()
            if line:
                contigs.append(line)
        except EOFError:
            break

    n50 = calculate_n_statistic(contigs, 50)
    n75 = calculate_n_statistic(contigs, 75)

    print(f"{n50} {n75}")

if __name__ == "__main__":
    main()

