## Usgae split_fasta_by_name.py <input.fasta> <names_to_remove.txt> > <filtered.fasta>

#!/usr/bin/env python3

from Bio import SeqIO
import sys

def main():
    if len(sys.argv) != 4:
        print("Usage: split_fasta_by_name.py <input.fasta> <names_to_remove.txt> > <filtered.fasta>", file=sys.stderr)
        sys.exit(1)
    
    input_fasta = sys.argv[1]
    names_to_remove_file = sys.argv[2]

    header_set = set(line.strip() for line in open(names_to_remove_file))

    with open(input_fasta, "r") as fasta_file:
        for seq_record in SeqIO.parse(fasta_file, "fasta"):
            try:
                header_set.remove(seq_record.name)
            except KeyError:
                print(seq_record.format("fasta"))
                continue

    if len(header_set) != 0:
        print(len(header_set), 'of the headers from the list were not identified in the input fasta file.', file=sys.stderr)

if __name__ == "__main__":
    main()
