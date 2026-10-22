import numpy as np
import pandas as pd
from scipy import stats

record_ks = {}
record_ka = {}
record_len = {}

# 1. Parse the human_chim_mouse_rat.txt file
with open('human_chim_mouse_rat.txt', 'r') as FILE:
    for line in FILE:
        line = line.strip()
        table = line.split('\t')
        try:
            # Based on the human_chim_mouse_rat.txt header:
            # Col 13 = Symbol, Col 24 = Ka (Human), Col 25 = Ks (Human), Col 15 = AA Compared (Length)
            gene_symbol = table[13]
            record_ks[gene_symbol] = float(table[25])  # ks value
            record_ka[gene_symbol] = float(table[24])  # ka value
            record_len[gene_symbol] = float(table[15]) # AA Compared (proxy for gene length)
        except:
            # Skips header or malformed lines
            pass

# 2. Parse the gene_consolidated_mask.txt file
data = []
with open('gene_consolidated_mask.txt', 'r') as FILE:
    for line in FILE:
        line = line.strip()
        table = line.split('\t')
        gene_symbol = table[0]
        
        if gene_symbol in record_ks:
            try:
                # Use table[2] for the 'predicted' column
                pred_val = float(table[2]) 
                data.append({
                    'Gene': gene_symbol,
                    'Ks': record_ks[gene_symbol],
                    'Ka': record_ka[gene_symbol],
                    'Length': record_len[gene_symbol],
                    'Value': pred_val
                })
            except:
                pass

# Convert to Pandas DataFrame for easy grouping
df = pd.DataFrame(data)

# 3. Create Length Groups
# pd.qcut divides the data into 4 quartiles so each group has roughly the same number of genes
df['Length_Group'] = pd.qcut(df['Length'], q=4, labels=['Short', 'Medium-Short', 'Medium-Long', 'Long'])

# 4. Calculate correlations per group
results = []
for group_name, group_df in df.groupby('Length_Group', observed=False):
    # We need at least 2 data points to calculate a correlation
    if len(group_df) < 2:
        continue
        
    ks_vect = group_df['Ks'].values
    ka_vect = group_df['Ka'].values
    pred_vect = group_df['Value'].values
    
    # Pearson correlations
    pearson_ks_pred = np.corrcoef(ks_vect, pred_vect)[0, 1]
    pearson_ka_pred = np.corrcoef(ka_vect, pred_vect)[0, 1]
    
    # Spearman correlations
    spearman_ks_pred, p_ks_pred = stats.spearmanr(ks_vect, pred_vect)
    spearman_ka_pred, p_ka_pred = stats.spearmanr(ka_vect, pred_vect)
    spearman_ka_ks, p_ka_ks = stats.spearmanr(ka_vect, ks_vect)
    
    results.append({
        'Length_Group': group_name,
        'N_Genes': len(group_df),
        'Min_Length': group_df['Length'].min(),
        'Max_Length': group_df['Length'].max(),
        'Pearson_Ks_vs_Pred': pearson_ks_pred,
        'Pearson_Ka_vs_Pred': pearson_ka_pred,
        'Spearman_Ks_vs_Pred': spearman_ks_pred,
        'Spearman_Ka_vs_Pred': spearman_ka_pred,
        'Spearman_Ka_vs_Ks': spearman_ka_ks
    })

# 5. Output and save the table
results_df = pd.DataFrame(results)

# Print the table to the console
print("\n--- Three Species Correlation Results Grouped by Gene Length (Predicted) ---")
print(results_df.to_string(index=False))

# Save the table to a file
output_filename = 'correlation_threespecies_by_length_groups_predicted.txt'
results_df.to_csv(output_filename, sep='\t', index=False)
print(f"\nResults successfully saved to '{output_filename}'")
