def longest_common_substring(strings):
    if not strings:
        return ""
    
    shortest_string = min(strings, key=len)
    longest_substring = ""
    
    for length in range(len(shortest_string), 0, -1):
        for start in range(len(shortest_string) - length + 1):
            substring = shortest_string[start:start+length]
            is_common = all(substring in string for string in strings)
            if is_common and len(substring) > len(longest_substring):
                longest_substring = substring
    
    return longest_substring

def main(input_file_path):
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
        
    longest_substring = longest_common_substring(strings)
    return longest_substring

# Example usage
input_file_path = "input.fasta"  # Replace with your input file path
longest_substring = main(input_file_path)
print(longest_substring)

