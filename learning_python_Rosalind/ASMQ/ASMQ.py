def calculate_n_statistic(contigs, percent):
    total_length = sum(len(contig) for contig in contigs)
    threshold_length = total_length * percent // 100

    contigs.sort(key=len, reverse=True)
    cumulative_length = 0
    for contig in contigs:
        cumulative_length += len(contig)
        if cumulative_length >= threshold_length:
            return len(contig)

def calculate_l_statistic(contigs, l_value):
    contigs.sort(key=len, reverse=True)
    cumulative_length = 0
    for idx, contig in enumerate(contigs, start=1):
        cumulative_length += len(contig)
        if cumulative_length >= l_value:
            return idx


def process_file(filename):
    contigs = []
    with open(filename, 'r') as file:
        for line in file:
            line = line.strip()
            if line:
                contigs.append(line)

    n50 = calculate_n_statistic(contigs, 50)
    n75 = calculate_n_statistic(contigs, 75)
    total_genome_size = sum(len(contig) for contig in contigs)
    l50 = calculate_l_statistic(contigs, total_genome_size * 0.5)

    return n50, n75, l50, total_genome_size

def main():
    input_filename = input("path to input file")
    output_filename = input("path to output file")

    with open(input_filename, 'r') as input_file, open(output_filename, 'w') as output_file:
        output_file.write("File Name\tN50\tN75\tL50\tTotal Genome Size\n")

        for line in input_file:
            dna_file = line.strip()
            n50, n75, l50, genome_size = process_file(dna_file)
            output_file.write(f"{dna_file}\t{n50}\t{n75}\t{l50}\t{genome_size}\n")

    print("Statistics written to", output_filename)

if __name__ == "__main__":
    main()
