import pandas as pd
import os

# Assuming the current directory contains the txt files
current_dir = os.getcwd()

# Read the txt files
observed_file = os.path.join(current_dir, "gene_observed_mask.txt")
predicted_file = os.path.join(current_dir, "gene_predicted_mask.txt")

# Read the files into pandas DataFrames
observed_df = pd.read_csv(observed_file, sep="\t", header=0)
predicted_df = pd.read_csv(predicted_file, sep="\t", header=0)

# Assuming the first column is the key for merging
key_column = observed_df.columns[0]

# Merge the dataframes based on the first column
merged_df = pd.merge(observed_df, predicted_df, on=key_column, suffixes=('_observed', '_predicted'))

# Compare the values for shared columns (excluding the key column)
shared_columns = [col for col in observed_df.columns if col in predicted_df.columns and col != key_column]

for col in shared_columns:
    merged_df[f'{col}_match'] = merged_df[f'{col}_observed'] == merged_df[f'{col}_predicted']

# Calculate the percentage of matching values for each shared column
match_percentages = merged_df[[f'{col}_match' for col in shared_columns]].mean() * 100

print("Percentage of matching values for shared columns:")
for col, percentage in match_percentages.items():
    original_col = col.replace('_match', '')
    print(f"{original_col}: {percentage:.2f}%")

# Save the merged and compared data to a new file
output_file = os.path.join(current_dir, "merged_and_compared.txt")
merged_df.to_csv(output_file, sep="\t", index=False)

print(f"\nMerged and compared data saved to: {output_file}")

