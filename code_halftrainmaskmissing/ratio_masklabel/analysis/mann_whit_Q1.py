from scipy.stats import mannwhitneyu
import sys
import numpy as np

# Function to extract lowest 25% (quartile 1) values
def get_lower_quartile(data):
    data_sorted = np.sort(data)  # Sort data
    lower_quartile = data_sorted[:len(data_sorted) // 4]  # Take the lowest 25%
    return lower_quartile

# Read data from file
def read_data(file_path):
    with open(file_path, 'r') as file:
        data = [float(line.strip()) for line in file]
    return data

# Read data from input files
file1_path = sys.argv[1]
file2_path = sys.argv[2]

group1 = read_data(file1_path)
group2 = read_data(file2_path)

# Extract lowest quartile (Q1) values
group1_q1 = get_lower_quartile(group1)
group2_q1 = get_lower_quartile(group2)

# Perform Mann-Whitney U test on lowest quartile values
stat, p_value = mannwhitneyu(group1_q1, group2_q1, alternative='two-sided')

# Print results
print(f"Mann-Whitney U statistic (Q1): {stat}")
print(f"P-value (Q1): {p_value}")

