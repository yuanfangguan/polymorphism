import pandas as pd
import random
import sys

# Step 1: Load the input and output file paths
input_file = sys.argv[1]  # Replace with your input file name
output_file = sys.argv[2]  # Replace with your output file name

# Load the file into a DataFrame
data = pd.read_csv(input_file, sep='\t', header=None, names=['chromosome', 'position', 'value', 'sequence'])

# Step 2: Function to shuffle the DNA sequence
def shuffle_sequence(sequence):
    sequence_list = list(sequence)
    random.shuffle(sequence_list)
    return ''.join(sequence_list)

# Step 3: Generate 100 shuffled sequences for each original sequence
simulated_data = []

for index, row in data.iterrows():
    original_sequence = row['sequence']
    shuffled_sequences = [shuffle_sequence(original_sequence) for _ in range(100)]
    # Append the chromosome, position, value, and tab-separated shuffled sequences
    for the_sequence in shuffled_sequences:
        simulated_data.append([row['chromosome'], row['position'], row['value'], the_sequence])


# Step 4: Save the simulated sequences to an output file
columns = ['chromosome', 'position', 'value','shuffled_seq']
simulated_df = pd.DataFrame(simulated_data, columns=columns)
simulated_df.to_csv(output_file, sep='\t', index=False)

print(f"Simulated sequences saved to {output_file}")

