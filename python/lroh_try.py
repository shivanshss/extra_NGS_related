import matplotlib.pyplot as plt

def read_lroh_file(lroh_file):
    lroh_data = []
    with open(lroh_file, 'r') as f:
        for line in f:
            # Skip header lines
            if line.startswith("#"):
                continue
            # Split the line into columns
            columns = line.strip().split()
            # Check if the line has the required number of columns
            if len(columns) < 3:
                continue
            # Extract start and end positions of the LROH segment
            try:
                start = int(columns[1])
                end = int(columns[2])
                lroh_data.append((start, end))
            except ValueError:
                # Skip lines with invalid data
                continue
    return lroh_data


def plot_lroh(lroh_data):
    fig, ax = plt.subplots(figsize=(10, 5))
    for start, end in lroh_data:
        ax.axvspan(start, end, alpha=0.3, color='blue')
    ax.set_xlabel("Position")
    ax.set_ylabel("LROH")
    ax.set_title("Long Runs of Homozygosity")
    plt.show()

if __name__ == "__main__":
    lroh_file = "/media/bloodmark/HDD6_SS_extra/w_newref/5.vcfstats/vcfstats_filter4/out.LROH"  # Replace with the path to your .LROH file
    lroh_data = read_lroh_file(lroh_file)
    plot_lroh(lroh_data)

