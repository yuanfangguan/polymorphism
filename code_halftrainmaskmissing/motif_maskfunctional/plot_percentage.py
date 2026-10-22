import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import gaussian_kde
plt.rcParams.update({'font.size': 14})

def calculate_at_percentage(sequence):
    count_a = sequence.count('A')
    count_t = sequence.count('T')
    total_bases = len(sequence.strip())
    at_percentage = (count_a + count_t) / total_bases * 100 if total_bases > 0 else 0
    return at_percentage

def process_file(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()
    
    at_percentages = []
    current_seq = ""
    for line in lines:
        line = line.strip()
        if line.startswith("NC_"):
            if current_seq:
                at_percentage = calculate_at_percentage(current_seq)
                at_percentages.append(at_percentage)
            chrom, pos, pred, current_seq = line.split('\t')
        else:
            current_seq += line
    
    if current_seq:
        at_percentage = calculate_at_percentage(current_seq)
        at_percentages.append(at_percentage)
    
    return at_percentages

# Paths to your two files
file_high = 'top_100000_windows.txt'
file_low = 'lowest_100000_windows.txt'

# Process both files
at_percentages_high = process_file(file_high)
at_percentages_low = process_file(file_low)

# Convert the percentages into numpy arrays for better handling
at_percentages_high = np.array(at_percentages_high)
at_percentages_low = np.array(at_percentages_low)

# Create density functions using gaussian_kde
density_high = gaussian_kde(at_percentages_high)
density_low = gaussian_kde(at_percentages_low)

# Generate a range of values for plotting the density curves
x_vals = np.linspace(0, 100, 500)  # Adjust the range based on your data

# Plot the density curves
plt.figure(figsize=(8, 5))

# Plot for the high file
plt.plot(x_vals, density_high(x_vals), label='Top 100,000 segments', color='blue')
plt.fill_between(x_vals, density_high(x_vals), color='blue', alpha=0.3)  # Fill with transparency

# Plot for the low file
plt.plot(x_vals, density_low(x_vals), label='Lowest 100,000 segments', color='red')
plt.fill_between(x_vals, density_low(x_vals), color='red', alpha=0.3)  # Fill with transparency

# Add labels, title, and legend
plt.xlabel('A/T Percentage')
plt.ylabel('Density')
plt.legend(loc='upper right')

# Show and save the plot
output_file = 'at_percentage_density_curve.png'
plt.savefig(output_file, dpi=300)
plt.grid(True)
plt.show()

