import numpy as np

# Load the .npy file
file_path = '../../../expression/NC_000001.11.result.npy'  # Replace with your file path
import numpy as np
import random

# Load the .npy file
data = np.load(file_path)

# Ensure the data has at least 100 * 128 elements to select from
if len(data) < 128 * 100:
        print("Data is too short to select 100 sections of length 128.")
else:
            # Randomly select 100 starting indices for the sections
            max_start_index = len(data) - 128  # Last valid start for a 128-length section
            start_indices = random.sample(range(max_start_index + 1), 100)

                        # Calculate the median for each section and print
            medians = []
            for start in start_indices:
                    section = data[start:start + 128]
                    section_median = np.median(section)
                    medians.append(section_median)
                    print(f"{section_median}")

data = np.load(file_path)

# Calculate the average
average = np.median(data)

print(f"The average of the .npy file is: {average}")

