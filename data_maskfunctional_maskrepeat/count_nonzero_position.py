import os
import pickle
import numpy as np

# Function to process one-hot encoded files and count masked and unmasked positions
def count_masked_and_unmasked(directory):
    total_masked = 0
    total_unmasked = 0

    # Iterate over all pkl files in the directory
    for filename in os.listdir(directory):
        if filename.endswith(".pkl"):
            with open(os.path.join(directory, filename), 'rb') as f:
                # Load the one-hot encoded numpy array
                one_hot_array = pickle.load(f)
                
                # Count the positions where all 4 elements are zeros (masked)
                masked_count = np.sum(np.all(one_hot_array == 0, axis=1))
                
                # Count the unmasked positions (those that are not all zeros)
                unmasked_count = one_hot_array.shape[0] - masked_count
                
                # Add to the totals
                total_masked += masked_count
                total_unmasked += unmasked_count

    return total_masked, total_unmasked

# Example usage
directory = './'
total_masked, total_unmasked = count_masked_and_unmasked(directory)

print(f'Total masked positions: {total_masked}')
print(f'Total unmasked positions: {total_unmasked}')

