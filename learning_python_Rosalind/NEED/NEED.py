#Given: Two GenBank IDs.

#Return: The maximum global alignment score between the DNA strings associated with these IDs. 

import subprocess

def get_alignment_score(genbank_id1, genbank_id2):
    # Fetch DNA sequences using GenBank API or other methods
    
    # Example DNA sequences (replace with actual sequences)
    sequence1 = "ACGT..."
    sequence2 = "TGCA..."

    # Run EMBOSS Needle tool
    command = [
        "needle",
        "-asequence", "stdin",
        "-bsequence", "stdin",
        "-gapopen", "10",
        "-gapextend", "1",
        "-datafile", "DNAfull",
        "-auto",
        "-stdout"
    ]

    process = subprocess.Popen(command, stdin=subprocess.PIPE, stdout=subprocess.PIPE, text=True)
    stdout, _ = process.communicate(input=f">{genbank_id1}\n{sequence1}\n>{genbank_id2}\n{sequence2}\n")

    # Parse alignment output to extract alignment score
    alignment_score = None
    for line in stdout.split('\n'):
        if "Alignment score" in line:
            alignment_score = float(line.split(":")[1].strip())
            break

    return alignment_score

genbank_id1 = "GenBankID1"
genbank_id2 = "GenBankID2"
max_alignment_score = get_alignment_score(genbank_id1, genbank_id2)
print("Maximum Alignment Score:", max_alignment_score)

