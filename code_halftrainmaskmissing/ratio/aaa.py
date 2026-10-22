import numpy as np
import random

def print_random_sections(npy_file, section_size=10, num_sections=5):
    # Step 1: Load the 1D npy file
    data = np.load(npy_file)
    
    # Step 2: Ensure the npy file is long enough for the sections
    data_length = len(data)
    if data_length < section_size:
        print(f"Data is too short. Length of data: {data_length}")
        return
    
    # Step 3: Randomly select start indices for the sections
    max_start_index = data_length - section_size
    random_start_indices = random.sample(range(max_start_index), num_sections)

    # Step 4: Print out the selected sections
    for i, start_idx in enumerate(random_start_indices):
        section = data[start_idx:start_idx + section_size]
        print(f"Section {i+1}: {section}")

# Usage
npy_file = 'processed_NC_000014.9.npy'  # Replace with your npy file path
print_random_sections(npy_file, section_size=10, num_sections=5)

