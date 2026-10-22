import os
import numpy as np
import matplotlib.pyplot as plt

# Directory containing the files - set this to the correct path
directory = './'
plt.rcParams.update({'font.size': 14})

# Define symbols and colors for each fold for differentiation
symbols = ['o', 's', 'D', 'v', '^']  # Different symbols for each fold
colors = ['blue', 'green', 'red', 'purple', 'orange']  # Different colors for each fold

# Define lists to store all train and test values
all_train_values = []
all_test_values = []
all_folds = []

# Loop through each fold (0 to 4) and each model within the fold (0 to 4)
for fold in range(5):  # There are 5 cross-validation folds
    for model in range(5):  # Each fold has 5 models

        # Construct the filenames
        train_file = f'../base_excludeN_5fold/cor.sep.{fold}.{model}'
        test_file = f'cor.sep.{fold}.{model}'

        # Construct the full file paths
        train_path = os.path.join(directory, train_file)
        test_path = os.path.join(directory, test_file)

        # Check if the files exist before loading
        if not os.path.isfile(train_path) or not os.path.isfile(test_path):
            print(f"File {train_path} or {test_path} not found. Skipping.")
            continue

        # Load the train and test correlation values
        train_correlations = np.loadtxt(train_path)
        test_correlations = np.loadtxt(test_path)

        # Extract the lower triangle of the correlation matrices (excluding the diagonal)
        train_corr_values = train_correlations[np.tril_indices_from(train_correlations, k=-1)]
        test_corr_values = test_correlations[np.tril_indices_from(test_correlations, k=-1)]

        # Append these values to the overall list
        all_train_values.extend(train_corr_values)
        all_test_values.extend(test_corr_values)
        all_folds.extend([fold] * len(train_corr_values))

# Convert to numpy arrays
all_train_values = np.array(all_train_values)
all_test_values = np.array(all_test_values)
all_folds = np.array(all_folds)

# Plotting the scatter plot for correlation values grouped by fold
plt.figure(figsize=(5.5, 5))
for fold in range(5):
    # Get the values for the current fold
    train_vals_fold = all_train_values[all_folds == fold]
    test_vals_fold = all_test_values[all_folds == fold]
    
    # Plot with a different symbol and color for each fold
    plt.scatter(train_vals_fold, test_vals_fold, alpha=0.6, 
                color=colors[fold], marker=symbols[fold], label=f'Fold {fold}')

plt.ylabel('Remove functional regions')
plt.xlabel('Base model')
plt.title('Correlation,individual model')

# Set the range [0-1] for both axes
plt.plot([0, 1], [0, 1], color='black', linestyle='--', linewidth=1, label='y = x')

plt.xlim(0.5, 0.7)
plt.ylim(0.5, 0.7)

# Remove the grid
plt.grid(False)
plt.tight_layout()

# Add legend
plt.legend()

# Save the plot
plt.savefig('correlation_scatter_plot_by_fold_compare_base.png', dpi=300)
plt.show()

