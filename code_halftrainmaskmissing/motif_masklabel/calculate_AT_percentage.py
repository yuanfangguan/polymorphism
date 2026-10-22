import sys
def calculate_at_percentage(sequence):

    # Count the number of A's and T's in the sequence
    count_a = sequence.count('A')
    count_t = sequence.count('T')
    
    # Calculate the total number of bases (excluding line breaks)
    total_bases = len(sequence.strip())
    
    # Calculate the A/T percentage
    at_percentage = (count_a + count_t) / total_bases * 100 if total_bases > 0 else 0
    return at_percentage

# Function to process the file and calculate A/T percentage for each line
def process_file(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()
    
    results = []
    current_seq = ""
    for line in lines:
        line = line.strip()
        if line.startswith("NC_"):
            # If this is a new record, process the previous one
            if current_seq:
                at_percentage = calculate_at_percentage(current_seq)
                results.append((chrom, pos, pred, at_percentage))
            # Start new record
            chrom, pos, pred, current_seq = line.split('\t')
        else:
            # Concatenate sequence if split across lines
            current_seq += line
    
    # Don't forget to process the last sequence
    if current_seq:
        at_percentage = calculate_at_percentage(current_seq)
        results.append((chrom, pos, pred, at_percentage))
    
    return results

# Path to your file
file_path = sys.argv[1]

# Process the file and print the results
results = process_file(file_path)
for chrom, pos, pred, at_percentage in results:
    print(f"{chrom}\t{pos}\t{pred}\t{at_percentage:.2f}% A/T content")

