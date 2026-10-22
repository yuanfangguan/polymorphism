import os
import numpy as np
from scipy.stats import mannwhitneyu

# Define directories
dir1 = 'predicted_output_segments'
dir2 = 'snp_segment_values'

# Output file
output_file = 'mann_whitney_results.txt'

with open(output_file, 'w') as out:
    out.write('Feature\tU_statistic\tp_value\n')
    
    # Loop over files in dir1
    for filename in os.listdir(dir1):
        file1_path = os.path.join(dir1, filename)
        file2_path = os.path.join(dir2, filename)
        
        # Check if corresponding file exists
        if not os.path.isfile(file2_path):
            print(f"Skipping {filename}, no corresponding file in {dir2}")
            continue
        
        # Load data
        data1 = np.loadtxt(file1_path)
        data2 = np.loadtxt(file2_path)
        
        # Perform Mann-Whitney U test
        try:
            stat, p = mannwhitneyu(data1, data2, alternative='two-sided')
            out.write(f"{filename}\t{stat}\t{p}\n")
        except Exception as e:
            print(f"Error processing {filename}: {e}")

print(f"Results saved to {output_file}")

