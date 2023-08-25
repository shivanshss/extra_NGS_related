from Bio import SeqIO
from Bio.Align import MultipleSeqAlignment
from Bio.Align import AlignInfo
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import os

def generate_combined_consensus(abi_files, output_dir):
    for row_num, (abi_file1, abi_file2) in enumerate(abi_files, start=1):
        # Parse the ABI files and extract sequences
        sequences = []
        for abi_file in [abi_file1, abi_file2]:
            with open(abi_file, "rb") as handle:
                record = SeqIO.read(handle, "abi")
                sequences.append(record.seq)

        # Create a multiple sequence alignment
        alignment = MultipleSeqAlignment([SeqRecord(seq) for seq in sequences])

        # Generate a consensus sequence
        summary_align = AlignInfo.SummaryInfo(alignment)
        consensus = summary_align.dumb_consensus()

        # Create a directory for each row if it doesn't exist
        row_dir = os.path.join(output_dir, f"row_{row_num}")
        os.makedirs(row_dir, exist_ok=True)

        # Save the consensus sequence in FASTA format
        output_file = os.path.join(row_dir, "consensus.fasta")
        with open(output_file, "w") as handle:
            record = SeqRecord(Seq(str(consensus)), id=f"row_{row_num}_consensus")
            SeqIO.write(record, handle, "fasta")

if __name__ == "__main__":
    input_filename = "input_list.txt"
    output_directory = "output_consensus"
    
    abi_files = []
    with open(input_filename, "r") as input_file:
        for line in input_file:
            file1, file2 = line.strip().split()
            abi_files.append((file1, file2))
    
    generate_combined_consensus(abi_files, output_directory)
