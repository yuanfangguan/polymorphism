import os
import numpy as np
import matplotlib.pyplot as plt

# Directory containing the files - set this to the correct path
directory = './'  # Adjust to the directory where your files are located
plt.rcParams.update({'font.size': 14})

# Define lists to store AUC values for training and test NNwise folds
train_auc_values_nn = []
test_auc_values_nn = []

# Loop through each fold to read AUC values
for i in range(5):  # Assuming there are 5 folds (0 to 4)
    # File paths
    train_auc_file_nn = f'../base_excludeN_5fold/auc.fold.NNwise{i}'
    test_auc_file_nn = f'auc.fold.NNwise{i}'

    # Construct the full file paths
    train_auc_path_nn = os.path.join(directory, train_auc_file_nn)
    test_auc_path_nn = os.path.join(directory, test_auc_file_nn)

    # Check if the AUC files exist before loading
    if os.path.isfile(train_auc_path_nn):
        train_auc_value = np.loadtxt(train_auc_path_nn)
        train_auc_values_nn.append(train_auc_value)
    else:
        print(f"Train AUC file {train_auc_path_nn} not found. Skipping.")

    if os.path.isfile(test_auc_path_nn):
        test_auc_value = np.loadtxt(test_auc_path_nn)
        test_auc_values_nn.append(test_auc_value)
    else:
        print(f"Test AUC file {test_auc_path_nn} not found. Skipping.")

# Convert lists to numpy arrays
train_auc_values_nn = np.array(train_auc_values_nn)
test_auc_values_nn = np.array(test_auc_values_nn)

# Plotting the scatter plot for AUC values
plt.figure(figsize=(5.5, 5))
plt.scatter(train_auc_values_nn, test_auc_values_nn, alpha=0.7, color='blue', marker='o', edgecolor='k',label='NN-wise AUC' )
plt.ylabel('Remove functional regions')
plt.xlabel('Base model')
plt.title('AUC at per nucleotide resolution')

#Test Set the range [0-1] for both axes
plt.plot([0, 1], [0, 1], color='black', linestyle='--', linewidth=1, label='y = x')

plt.xlim(0.70, 0.73)
plt.ylim(0.70, 0.73)
plt.tight_layout()

# Remove the grid
plt.grid(False)

# Add legend
plt.legend()

# Save the plot
plt.savefig('auc_comparison_compare_base.png', dpi=300)
plt.show()

