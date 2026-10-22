import numpy as np
from scipy.stats import ttest_ind

import sys
# Load the data from the two files
file1 = sys.argv[1]  # Replace with your actual file path
file2 = sys.argv[2]  # Replace with your actual file path

# Read the data from each file
data1 = np.loadtxt(file1)
data2 = np.loadtxt(file2)

# Perform the t-test (independent two-sample t-test)
t_stat, p_value = ttest_ind(data1, data2)

# Print the results
print(f"T-statistic: {t_stat}")
print(f"P-value: {p_value}")

