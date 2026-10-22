import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# 1. Set plot style for publication (Increased global font_scale from 1.2 to 1.5)
sns.set_theme(style="whitegrid", context="paper", font_scale=1.5)
plt.rcParams['font.family'] = 'sans-serif'

# 2. File paths (update these if your file names are different)
files = {
    'pred_chimp': 'correlation_by_length_groups_predicted_chimp_ks.txt',
    'obs_chimp': 'correlation_by_length_groups_chimp_observed_ks.txt',
    'pred_three': 'correlation_threespecies_by_length_groups_predicted.txt',
    'obs_three': 'correlation_threespecies_by_length_groups_observed.txt'
}

# 3. Load data
dfs = {name: pd.read_csv(filepath, sep='\t') for name, filepath in files.items()}

# Ensure the Length_Group is treated as a categorical variable in the correct order
length_order = ['Short', 'Medium-Short', 'Medium-Long', 'Long']
for df in dfs.values():
    df['Length_Group'] = pd.Categorical(df['Length_Group'], categories=length_order, ordered=True)

# 4. Create the 2x2 figure (Increased figsize from 12,10 to 14,12 to accommodate larger fonts)
fig, axes = plt.subplots(2, 2, figsize=(14, 12))

def plot_panel(ax, df, is_predicted, title):
    # Determine column names based on whether it's predicted or observed
    ks_col = 'Spearman_Ks_vs_Pred' if is_predicted else 'Spearman_Ks_vs_Obs'
    ka_col = 'Spearman_Ka_vs_Pred' if is_predicted else 'Spearman_Ka_vs_Obs'
    
    # Plot D_s (formerly Ks) - slightly increased markersize/linewidth for balance
    ax.plot(df['Length_Group'], df[ks_col], 
            marker='o', markersize=9, linewidth=3.0, 
            color='#1f77b4', label=r'vs. $D_s$ (Neutral)')
    
    # Plot D_n (formerly Ka) - slightly increased markersize/linewidth for balance
    ax.plot(df['Length_Group'], df[ka_col], 
            marker='s', markersize=9, linewidth=3.0, 
            color='#d62728', label=r'vs. $D_n$ (Selection)')
    
    # Formatting (Increased font sizes)
    ax.set_title(title, fontsize=18, fontweight='bold', pad=14)
    ax.set_ylabel("Spearman Correlation", fontsize=16)
    ax.set_xlabel("Gene Length Quartiles", fontsize=16)
    
    # Fix Y-axis to be the same across all plots for fair visual comparison
    ax.set_ylim(-0.15, 0.50) 
    
    # Add a horizontal line at 0 for reference
    ax.axhline(0, color='black', linestyle='--', linewidth=1.5, alpha=0.7)
    
    # Increased legend font size
    ax.legend(loc='upper left', frameon=True, fontsize=14)
    
    # Increased tick label sizes for both X and Y axes
    ax.tick_params(axis='x', labelsize=14, labelrotation=15)
    ax.tick_params(axis='y', labelsize=14)

# 5. Plot each panel 
# Top Row: Chimpanzee
plot_panel(axes[0, 0], dfs['pred_chimp'], is_predicted=True, title="E. Predicted Mutation (Chimp Outgroup)")
plot_panel(axes[0, 1], dfs['obs_chimp'], is_predicted=False, title="F. Observed Mutation (Chimp Outgroup)")

# Bottom Row: Mouse/Rat
plot_panel(axes[1, 0], dfs['pred_three'], is_predicted=True, title="G. Predicted Mutation (Mouse/Rat Outgroup)")
plot_panel(axes[1, 1], dfs['obs_three'], is_predicted=False, title="H. Observed Mutation (Mouse/Rat Outgroup)")

# 6. Final layout adjustments and save
plt.tight_layout()
output_fig = "Correlation_Length_Trends.pdf"
plt.savefig(output_fig, dpi=300, bbox_inches='tight')
print(f"Figure successfully saved as {output_fig}")

# Show the plot if running in an interactive environment (like Jupyter)
plt.show()
