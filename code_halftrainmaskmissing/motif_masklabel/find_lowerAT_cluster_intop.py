def calculate_at_percentage(sequence):
    """
    Calculate the percentage of A/T content in a given DNA sequence.
    """
    count_a = sequence.count('A')
    count_t = sequence.count('T')
    
    # Calculate the total number of bases (excluding line breaks)
    total_bases = len(sequence.strip())
    
    # Calculate the A/T percentage
    at_percentage = (count_a + count_t) / total_bases * 100 if total_bases > 0 else 0
    return at_percentage

def filter_low_at_content(input_file, output_file, threshold=30):
    """
    Read the input file, calculate A/T content, and write lines with A/T content below the threshold.
    """
    with open(input_file, 'r') as file:
        lines = file.readlines()
    
    with open(output_file, 'w') as outfile:
        for line in lines:
            parts = line.strip().split('\t')
            if len(parts) == 4:
                sequence = parts[3]
                at_percentage = calculate_at_percentage(sequence)
                
                # If A/T percentage is below the threshold, write to output file
                if at_percentage < threshold:
                    outfile.write(line)

# Specify input and output file paths
input_file = 'top_100000_windows.txt'  # Replace with the actual file path
output_file = 'top_100000_windows.txt_low_at_content.txt'

# Call the function to filter lines with AT content below 30%
filter_low_at_content(input_file, output_file)

