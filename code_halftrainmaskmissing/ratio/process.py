import numpy as np
import os
import pickle

def load_npy_file(file_path):
    """Load a single .npy file."""
    return np.load(file_path)

def load_pkl_file(file_path):
    """Load a single .pkl file (numerator)."""
    with open(file_path, 'rb') as f:
        return pickle.load(f)


def moving_window_ratio_single_pair(numerator_file, denominator_file, output_file, window_size=2048, step_size=128):
    """Compute the moving window ratio for a single pair of numerator .pkl and denominator .npy files."""
    
    # Load the numerator from .pkl and denominator from .npy
    numerator_data = load_pkl_file(numerator_file)
    denominator_data = load_npy_file(denominator_file)
    
    # Use the shortest length of the two arrays if they differ
    min_length = min(len(numerator_data), len(denominator_data))
    numerator_data = numerator_data[:min_length]
    denominator_data = denominator_data[:min_length]
    
    # Initialize the output array with NaN for labeling and clearer missing value handling
    output = np.full(min_length, np.nan)
    
    # Iterate over the data using a sliding window
    for start in range(0, len(numerator_data) - window_size + 1, step_size):
        end = start + window_size
        
        # Compute the sum for both windows
        window_numerator = numerator_data[start:end]
        window_denominator = denominator_data[start:end]
        
        # Handle cases where the sum in either window is 0
        with np.errstate(divide='ignore', invalid='ignore'):  # Ignore divide by zero warnings
            ratio = np.divide(window_numerator, window_denominator, where=window_denominator != 0)
            ratio[window_denominator == 0] = np.nan  # Assign NaN where denominator is zero
           # print('observed:',window_numerator[0:10])
           # print('predicted:',window_denominator[0:10])
           # print(ratio[0:10])

        # Assign the ratio to the corresponding positions in the output array
        output[start:end] = ratio
    
    # Save the result to a .npy file
    np.save(output_file, output)
    print(f"Processed and saved: {output_file}")
    print(len(output))




def process_folders_pairwise(numerator_folder, denominator_folder, output_folder, window_size=2048, step_size=128):
    """Process files in two folders one pair at a time and save the output."""
    
    # Get lists of files in both folders
    denominator_files = [f for f in os.listdir(denominator_folder) if f.endswith('.pred.npy')]
    numerator_files = [f for f in os.listdir(numerator_folder) if f.endswith('_label.pkl')]
    
    # Process each pair based on the matching file prefix
    for denom_file in denominator_files:
        # Extract the prefix (e.g., 'NC_000009.12') from the denominator file
        prefix = denom_file.replace('.pred.npy', '')
        
        # Find the corresponding numerator file (e.g., 'NC_000009.12_label.pkl')
        numerator_file = os.path.join(numerator_folder, f'{prefix}_label.pkl')
        if not os.path.exists(numerator_file):
            print(f"No matching numerator file found for {denom_file}")
            continue
        
        # Define paths for input files and output file
        denom_file_path = os.path.join(denominator_folder, denom_file)
        output_file_path = os.path.join(output_folder, f'processed_{prefix}.npy')
        
        # Process the current pair and save the result
        moving_window_ratio_single_pair(numerator_file, denom_file_path, output_file_path, window_size, step_size)

# Example usage
numerator_folder = '../../data'  # Folder containing numerator .pkl files
denominator_folder = '../../calibrated_predictions'  # Folder containing denominator .npy files
output_folder = './'  # Folder to save the results

# Process and save each file pair one at a time
process_folders_pairwise(numerator_folder, denominator_folder, output_folder)

