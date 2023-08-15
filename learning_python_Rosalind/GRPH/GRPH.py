def generate_overlap_graph(strings, k):
    adjacency_list = []
    
    for i in range(len(strings)):
        for j in range(len(strings)):
            if i != j and strings[i][-k:] == strings[j][:k]:
                adjacency_list.append((strings[i], strings[j]))
    
    return adjacency_list

def main(input_file_path, k):
    strings = []
    current_string = ""
    
    with open(input_file_path, 'r') as input_file:
        for line in input_file:
            line = line.strip()
            if line.startswith('>'):
                if current_string:
                    strings.append(current_string)
                current_string = ""
            else:
                current_string += line
        strings.append(current_string)
        
    adjacency_list = generate_overlap_graph(strings, k)
    
    return adjacency_list

# Example usage
input_file_path = "input.fasta"  # Replace with your input file path
k = 3
adjacency_list = main(input_file_path, k)

# Print the adjacency list
for edge in adjacency_list:
    print(edge[0], edge[1])

