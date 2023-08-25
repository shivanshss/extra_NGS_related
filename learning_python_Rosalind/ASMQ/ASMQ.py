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


def process_fasta(filename):
    contig_lengths = []
    current_sequence = ""

    with open(filename, 'r') as file:
        for line in file:
            line = line.strip()
            if line.startswith(">"):
                if current_sequence:
                    contig_lengths.append(len(current_sequence))
                    current_sequence = ""
            else:
                current_sequence += line

        if current_sequence:  # To handle the last sequence in the file
            contig_lengths.append(len(current_sequence))

    n50 = calculate_n_statistic(contig_lengths, 50)
    n75 = calculate_n_statistic(contig_lengths, 75)
    total_genome_size = sum(contig_lengths)
    l50 = calculate_l_statistic(contig_lengths, total_genome_size * 0.5)

    return n50, n75, l50, total_genome_size

def main():
    input_filename = input("Enter the input filename: ")
    output_filename = input("Enter the output filename: ")

    with open(input_filename, 'r') as input_file, open(output_filename, 'w') as output_file:
        output_file.write("File Name\tN50\tN75\tL50\tTotal Genome Size\n")

        for line in input_file:
            fasta_file = line.strip()
            n50, n75, l50, genome_size = process_fasta(fasta_file)
            output_file.write(f"{fasta_file}\t{n50}\t{n75}\t{l50}\t{genome_size}\n")

    print("Statistics written to", output_filename)

if __name__ == "__main__":
    main()
