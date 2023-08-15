def find_substring_locations(s, t):
    locations = []
    t_length = len(t)
    
    for i in range(len(s) - t_length + 1):
        if s[i:i+t_length] == t:
            locations.append(i + 1)  # Adding 1 to match 1-based indexing
        
    return locations

def process_files(s_file_path, t_file_path):
    with open(s_file_path, 'r') as s_file, open(t_file_path, 'r') as t_file:
        s_strings = s_file.read().splitlines()
        t_strings = t_file.read().splitlines()
        
    for s, t in zip(s_strings, t_strings):
        locations = find_substring_locations(s, t)
        print(" ".join(map(str, locations)))

# Example usage
s_file_path = "s_strings.txt"
t_file_path = "t_strings.txt"
process_files(s_file_path, t_file_path)

