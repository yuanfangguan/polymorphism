import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D

# Load the data
df = pd.read_csv('combined_comparison.txt', sep='\t')
print(df.columns)

# Sort the DataFrame based on the Observed values
df_sorted = df.sort_values(by='Observed')

# Split the data into two halves
mid_index = len(df_sorted) // 2
df_first_half = df_sorted[:mid_index]
df_second_half = df_sorted[mid_index:]

# Combine the two halves back together, ensuring there's no duplicate index
df_combined = pd.concat([df_first_half, df_second_half]).reset_index(drop=True)

# Add a column for the ratio of Observed to Predicted
df_combined['Ratio'] = df_combined['Observed'] / df_combined['Predicted']

# Classify the ratio into categories for coloring
df_combined['Color'] = df_combined['Ratio'].apply(
    lambda x: 'red' if x < 0.9 else ('blue' if 0.9 <= x <= 1.1 else 'green')
)

# Plot the combined data with an additional panel
def plot_combined_with_ratio(df, title, filename):
    fig = plt.figure(figsize=(20, 10))
    gs = GridSpec(2, 1, height_ratios=[1, 3])

    # Panel 1: Ratio (OVEP) with colored dots
    ax1 = fig.add_subplot(gs[0])
    index = range(len(df))
    for i, color in enumerate(df['Color']):
        ax1.plot(i, df['Ratio'].iloc[i], 'o', color=color)
    
    # Add a horizontal dashed line at y=1
    ax1.axhline(y=1, color='black', linestyle='--', linewidth=1)
    
    ax1.set_ylabel('OVEP',fontsize=14)# = #observed/#predicted', fontsize=14)
    ax1.set_xticks([])
    ax1.tick_params(axis='y', labelsize=12)

    legend_elements = [
        Line2D([0], [0], marker='o', color='w', markerfacecolor='red', markersize=10, label='OVEP < 0.9'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor='blue', markersize=10, label='0.9 <= OVEP <= 1.1'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor='green', markersize=10, label='OVEP > 1.1')
    ]

    ax1.legend(handles=legend_elements, fontsize=12, loc='upper left')
    #ax1.legend(['OVEP < 0.9 (red)', '0.9 <= OVEP <= 1.1 (blue)', 'OVEP > 1.1 (green)'], 
               #fontsize=12, loc='upper left')
    ax1.set_title('OVEP (observed relative to expected number of polymorphic sites) ', fontsize=14)

    # Panel 2: Bar plot
    ax2 = fig.add_subplot(gs[1])
    bar_width = 0.35
    opacity = 0.8

    rects1 = ax2.bar(index, df['Observed'], bar_width, alpha=opacity, color='red', label='Observed')
    rects2 = ax2.bar([p + bar_width for p in index], df['Predicted'], bar_width, alpha=opacity, color='blue', label='Predicted')

    df['Type'] = df['Type'].str.replace('_', ' ')
    ax2.set_xticks([p + bar_width / 2 for p in index])
    ax2.set_xticklabels(df['Type'], rotation=45, ha='right', fontsize=12)
    ax2.set_ylabel('Fraction of polymorphic sites\n', fontsize=14)
    ax2.legend(fontsize=12, loc='upper left')
    ax2.tick_params(axis='y', labelsize=12)

    max_value = max(df['Observed'].max(), df['Predicted'].max())
    ax2.set_ylim(0, max_value * 1.1)

    plt.tight_layout()
    plt.savefig(f'{filename}.png', dpi=300)
    plt.show()

# Plot and save the combined data with the additional panel
plot_combined_with_ratio(df_combined, '', 'combined_comparison_with_ratio')

