import numpy as np

# Path to the directory containing the .npy files
npy_dir = '../../../expression'

# Path to the text file containing the segments
segment_file = 'random_100_32_windows.txt'


# Function to calculate the average expression level for a segment
def calculate_average_expression(chromosome, start_position, npy_dir, segment_length=12800):
    # Load the npy file for the chromosome
    npy_file = f"{npy_dir}/{chromosome}.result.npy"
    try:
        coverage_array = np.load(npy_file)
    except FileNotFoundError:
        print(f"File not found for {chromosome}")
        return None

    # Extract the segment from the coverage array
    end_position = start_position + segment_length
    if end_position > len(coverage_array):
        print(f"End position exceeds the length of the chromosome for {chromosome}. Skipping.")
        return None
    
    segment_coverage = coverage_array[start_position:end_position]
    print(start_position,end_position,segment_coverage)

    # Calculate and return the average expression level
    avg_expression = np.median(segment_coverage)
    return avg_expression

# Read the segment file and process each line
with open(segment_file, 'r') as file:
    lines = file.readlines()

    for line in lines:
        columns = line.strip().split('\t')
        
        # Extract chromosome, start position, and sequence
        chromosome = columns[0]  # Chromosome name (e.g., "NC_000010.11")
        start_position = int(columns[1])  # Starting position (e.g., 93965376)
        
        # Calculate the average expression level for the segment
        avg_expression = calculate_average_expression(chromosome, start_position, npy_dir)
        
        if avg_expression is not None:
            print(f"Chromosome: {chromosome}, Start: {start_position}, Average Expression: {avg_expression}")


