import pandas as pd
import random

import sys
# Step 1: Load the input file
input_file = sys.argv[1]# Replace with your file name
output_file =sys.argv[2] 

# Load the file into a DataFrame
data = pd.read_csv(input_file, sep='\t', header=None, names=['chromosome', 'position', 'value', 'sequence'])

# Step 2: Function to shuffle the DNA sequence
def shuffle_sequence(sequence):
    sequence_list = list(sequence)
    random.shuffle(sequence_list)
    return ''.join(sequence_list)

# Step 3: Generate 10 shuffled sequences for each original sequence
simulated_data = []

for index, row in data.iterrows():
    original_sequence = row['sequence']
    for i in range(10):  # Generate 10 shuffled sequences
        shuffled_sequence = shuffle_sequence(original_sequence)
        simulated_data.append([row['chromosome'], row['position'], row['value'], shuffled_sequence])

# Step 4: Save the simulated sequences to an output file
simulated_df = pd.DataFrame(simulated_data, columns=['chromosome', 'position', 'value', 'simulated_sequence'])
simulated_df.to_csv(output_file, sep='\t', index=False)

print(f"Simulated sequences saved to {output_file}")

