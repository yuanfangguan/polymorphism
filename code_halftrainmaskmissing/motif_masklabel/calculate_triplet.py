from collections import defaultdict
import sys

def calculate_triplet_distribution(file_path):
    # Dictionary to store triplet counts
    triplet_counts = defaultdict(int)
    total_triplets = 0

    with open(file_path, 'r') as file:
        for line in file:
            # Skip lines that do not contain sequences
            columns = line.strip().split('\t')
            if len(columns) < 4:
                continue

            # Extract the sequence
            sequence = ''.join(columns[3:])  # Combine if split into multiple lines
            sequence = sequence.replace("\n", "").upper()  # Remove newlines and normalize case

            # Calculate triplet frequencies
            for i in range(len(sequence) - 2):
                triplet = sequence[i:i+3]
                if len(triplet) == 3:  # Ensure triplet length
                    triplet_counts[triplet] += 1
                    total_triplets += 1

    # Calculate frequencies
    triplet_frequencies = {triplet: count / total_triplets for triplet, count in triplet_counts.items()}

    return triplet_frequencies

# Input file path
file_path = sys.argv[1]  # Replace with your file's path

# Calculate and print the triplet distribution
triplet_distribution = calculate_triplet_distribution(file_path)
for triplet, frequency in sorted(triplet_distribution.items()):
    print(f"{triplet}: {frequency:.4f}")

