import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import mannwhitneyu

# Load the first file (single-column)
file1 = "64_mean_expression.txt"
data1 = []
with open(file1, "r") as f:
    for line in f:
        value = line.strip()
        data1.append(float(value) if value != "nan" else 0)  # Convert 'nan' to 0

# Load the second file (extract last column, treating missing values as 0)
file2 = "whole_genome/64average_expression_results.txt"
data2 = []
with open(file2, "r") as f:
    next(f)  # Skip header
    for line in f:
        columns = line.strip().split("\t")
        value = columns[-1] if len(columns) > 3 else "0"  # If last column is missing, set to 0
        data2.append(float(value) if value != "nan" else 0)

# Convert to numpy arrays and apply log(x+1) transformation
data1 = np.log1p(np.array(data1))  # log(x+1) transformation
data2 = np.log1p(np.array(data2))

# Mann-Whitney U test
u_statistic, p_value_top_quarter = mannwhitneyu(data1, data2)
print(f"P-value for the Mann-Whitney U test: {p_value_top_quarter:.4f}")

# Create a single density plot
plt.figure(figsize=(8, 5))
sns.kdeplot(data1, label="random sections (log)", color="blue", fill=True, alpha=0.5)
sns.kdeplot(data2, label="top predictions (log)", color="orange", fill=True, alpha=0.5)
plt.ylabel("Density")
plt.xlabel("Log(Avg Expression + 1)")
plt.title("Expression Distribution 64bp (Density)")
plt.legend()
plt.grid(True)

# Save and show the plot
plt.savefig("expression_density_comparison.png", dpi=300)
plt.show()

