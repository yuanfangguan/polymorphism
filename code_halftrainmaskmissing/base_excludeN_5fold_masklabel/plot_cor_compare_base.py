import os
import numpy as np
import matplotlib.pyplot as plt

# Directory containing the files - set this to the correct path
directory = './'
plt.rcParams.update({'font.size': 14})

# Define lists to store the correlation values for train and test
train_values = []
test_values = []
train_values_nn = []
test_values_nn = []

# Loop through each train.cor.fold file (excluding NNwise files)
for i in range(5):  # Assuming there are 5 folds (0 to 4)
    # Standard fold files
    train_file = f'../base_excludeN_5fold/cor.fold.{i}'
    test_file = f'cor.fold.{i}'

    # NNwise fold files
    train_file_nn = f'../base_excludeN_5fold/cor.fold.NNwise.{i}'
    test_file_nn = f'cor.fold.NNwise.{i}'

    # Construct the full file paths
    train_path = os.path.join(directory, train_file)
    test_path = os.path.join(directory, test_file)
    train_path_nn = os.path.join(directory, train_file_nn)
    test_path_nn = os.path.join(directory, test_file_nn)

    # Check if files exist before loading for standard folds
    if not os.path.isfile(train_path) or not os.path.isfile(test_path):
        print(f"File {train_path} or {test_path} not found. Skipping.")
        continue

    # Load the train and test correlation values for standard folds
    train_correlations = np.loadtxt(train_path)
    test_correlations = np.loadtxt(test_path)

    # Extract the lower triangle of the correlation matrices (excluding the diagonal)
    train_corr_values = train_correlations[np.tril_indices_from(train_correlations, k=-1)]
    test_corr_values = test_correlations[np.tril_indices_from(test_correlations, k=-1)]

    # Append these values to the overall list
    train_values.extend(train_corr_values)
    test_values.extend(test_corr_values)

    # Check if files exist before loading for NNwise folds
    if not os.path.isfile(train_path_nn) or not os.path.isfile(test_path_nn):
        print(f"File {train_path_nn} or {test_path_nn} not found. Skipping.")
        continue

    # Load the train and test correlation values for NNwise folds
    train_correlations_nn = np.loadtxt(train_path_nn)
    test_correlations_nn = np.loadtxt(test_path_nn)

    # Extract the lower triangle of the correlation matrices (excluding the diagonal)
    train_corr_values_nn = train_correlations_nn[np.tril_indices_from(train_correlations_nn, k=-1)]
    test_corr_values_nn = test_correlations_nn[np.tril_indices_from(test_correlations_nn, k=-1)]

    # Append these values to the overall list for NNwise
    train_values_nn.extend(train_corr_values_nn)
    test_values_nn.extend(test_corr_values_nn)

# Convert to numpy arrays
train_values = np.array(train_values)
test_values = np.array(test_values)
train_values_nn = np.array(train_values_nn)
test_values_nn = np.array(test_values_nn)

# Plotting the scatter plot
plt.figure(figsize=(5.5, 5))
plt.scatter(train_values, test_values, alpha=0.5, color='red', label='2048 Chunks')
plt.scatter(train_values_nn, test_values_nn, alpha=0.5, color='blue', marker='^', label='Nucleotide-wise')
plt.xlabel('Base model')
plt.ylabel('Remove labeled regions')
plt.title('Correlation coefficient')

# Set the range [0-1] for both axes
plt.plot([0, 1], [0, 1], color='black', linestyle='--', linewidth=1, label='y = x')

plt.xlim(0, 1)
plt.ylim(0, 1)

# Remove the grid
plt.grid(False)

# Add legend
plt.legend()

# Save the plot
plt.savefig('correlation_scatter_plot_compare_base.png', dpi=300)
plt.show()

