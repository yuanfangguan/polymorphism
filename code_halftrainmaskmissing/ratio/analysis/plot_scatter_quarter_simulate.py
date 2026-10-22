import os
import numpy as np
import matplotlib.pyplot as plt
import glob

# Directory containing the files - set this to the correct path if necessary
directory = './'

# Function to calculate the median of values in a file
def calculate_median(file_path):
    data = np.loadtxt(file_path)
    return np.median(data)

# Function to calculate the median of the lowest 25% values
def calculate_lower_quartile_median(file_path):
    data = np.loadtxt(file_path)
    data_sorted = np.sort(data)  # Sort the data
    lower_quartile = data_sorted[:len(data_sorted) // 4]  # Take the lowest 25%
    return np.median(lower_quartile)

# Define sizes we're interested in
sizes = [64, 128, 256, 512, 1024, 2048]

# Define markers for each size
markers = ['o', 's', 'D', '^', 'v', 'P']

# Lists to store the median values
random_medians = []
lowest_medians = []
random_q1_medians = []
lowest_q1_medians = []
size_labels = []

# Loop through each size and collect the median values
for size in sizes:
    # Construct the filenames
    random_file = f'100times_simulated_{size}.txt.viennafold.summary'
    lowest_file = f'lowest_100_{size}_windows.txt.viennafold.summary'

    # Construct the full file paths
    random_path = os.path.join(directory, random_file)
    lowest_path = os.path.join(directory, lowest_file)

    # Calculate medians if files exist
    if os.path.isfile(random_path) and os.path.isfile(lowest_path):
        random_median = calculate_median(random_path) / size
        lowest_median = calculate_median(lowest_path) / size

        random_q1_median = calculate_lower_quartile_median(random_path) / size
        lowest_q1_median = calculate_lower_quartile_median(lowest_path) / size

        random_medians.append(random_median)
        lowest_medians.append(lowest_median)
        random_q1_medians.append(random_q1_median)
        lowest_q1_medians.append(lowest_q1_median)
        size_labels.append(size)
    else:
        print(f"Files for size {size} not found. Skipping this size.")

# Convert lists to numpy arrays
random_medians = np.array(random_medians)
lowest_medians = np.array(lowest_medians)
random_q1_medians = np.array(random_q1_medians)
lowest_q1_medians = np.array(lowest_q1_medians)

# Plotting the scatter plot
plt.figure(figsize=(5.5, 5))

# Original medians (main dataset)
for i, size in enumerate(sizes):
    if size in size_labels:
        idx = size_labels.index(size)
        plt.scatter(random_medians[idx], lowest_medians[idx], alpha=0.8, 
                    color='blue', marker=markers[i], edgecolor='black', 
                    label=f'Size {size}')

        # Annotate the point with the size
        plt.text(random_medians[idx], lowest_medians[idx], f'{size}', fontsize=9,
                 verticalalignment='bottom', horizontalalignment='right')

# New series: Lowest 25% median values
for i, size in enumerate(sizes):
    if size in size_labels:
        idx = size_labels.index(size)
        plt.scatter(random_q1_medians[idx], lowest_q1_medians[idx], alpha=0.8, 
                    color='green', marker=markers[i], edgecolor='black', 
                    label=f'Q1 Size {size}')# if i == 0 else "")

        # Annotate the point with the size
        plt.text(random_q1_medians[idx], lowest_q1_medians[idx], f'{size}', fontsize=9,
                 verticalalignment='top', horizontalalignment='left')

# Add the y = x reference line
plt.plot([min(random_medians) - 0.3, max(random_medians) + 0.1], 
         [min(random_medians) - 0.3, max(random_medians) + 0.1], 
         linestyle='--', color='red', label='y = x')

plt.xlabel('Median free energy per nucleotide (Scrambled sequence)')
plt.ylabel('Median free energy per nucleotide (Predicted structure)')
plt.xlim(-0.4, -0.15)
plt.ylim(-0.4, -0.15)
plt.legend()
plt.tight_layout()

# Save the plot
plt.savefig('median_scatter_plot_with_q1_simulate.png', dpi=300)
plt.show()

