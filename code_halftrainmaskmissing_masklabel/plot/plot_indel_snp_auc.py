#!/usr/bin/env python3
import os
import glob
import pandas as pd
import matplotlib.pyplot as plt

# Define file paths and labels
file_paths = {
    "Train: all; Test: all": "../base_excludeN_5fold_masklabel/auc.fold.*",
    "Train: INDEL; Test: INDEL": "../base_excludeN_5fold_indelonly_masklabel/auc.fold.*",
    "Train: all; Test: INDEL": "../base_excludeN_5fold_indelonlyintest_masklabel/auc.fold.*",
    "Train: INDEL; Test: all": "../base_excludeN_5fold_indelonlyonlyintrain_masklabel/auc.fold.*",
    "Train: INDEL; Test: SNP": "../base_excludeN_5fold_indelonlyonlyintrain_nonindelonlyintest_masklabel/auc.fold.*",
    "Train: SNP; Test: SNP": "../base_excludeN_5fold_nonindel_masklabel/auc.fold.*",
    "Train: all; Test: SNP": "../base_excludeN_5fold_nonindelonlyintest_masklabel/auc.fold.*",
    "Train: SNP; Test: all": "../base_excludeN_5fold_nonindelonlyintrain_masklabel/auc.fold.*",
    "Train: SNP; Test: INDEL": "../base_excludeN_5fold_nonindelonlyintrain_indelonlyintest_masklabel/auc.fold.*"
}

# Function to load data excluding NNwise files
def load_data(files_pattern):
    all_data = []
    for file in sorted(glob.glob(files_pattern)):
#        if "NNwise" not in file:
            data = pd.read_csv(file, delim_whitespace=True, header=None)
            all_data.append(data.iloc[0, 0])
    return all_data

# Load data for each category
data_dict = {label: load_data(pattern) for label, pattern in file_paths.items()}

# Create a box plot
plt.figure(figsize=(9, 4.5))
plt.boxplot([data_dict[label] for label in file_paths.keys()], labels=file_paths.keys(), patch_artist=True)
plt.xticks(rotation=45, ha='right')
plt.ylabel("AUC", fontsize=14)
plt.xlabel("Data Sets", fontsize=14)
plt.title("Train/Test Configurations", fontsize=16)
plt.ylim(0.5,1)
plt.tick_params(axis='both', which='major', labelsize=14)
plt.tight_layout()

# Save the figure
plt.savefig("box_plot_indel_snp_configurations_auc_masklabel.png")
plt.show()

