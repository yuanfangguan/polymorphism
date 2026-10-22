import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Load the data
df = pd.read_csv('combined_comparison.txt', sep='\t')

# Filter rows where Observed / Predicted < 0.9
df_filtered = df[df['Observed'] / df['Predicted'] > 1.1]

# Define a list of marker styles and colors for the scatter plot
markers = ['o', 's', 'D', '^', 'v', '>', '<', 'p', '*', 'h', 'H', '+', 'x', '|']
colors = plt.cm.get_cmap('tab20', len(df_filtered['Type'].unique()))  # Generate a color map with a unique color for each type

# Sort the filtered dataframe by 'Observed'
df_filtered = df_filtered.sort_values(by='Observed').reset_index(drop=True)

# Create a scatter plot with different markers and a legend
def plot_scatter_with_legend(df, title, filename):
    fig, ax = plt.subplots(figsize=(9, 8))

    # Set axis labels and title with font size 14
    ax.set_xlabel('Predicted', fontsize=14)
    ax.set_ylabel('Observed', fontsize=14)
    ax.set_title(title, fontsize=14)

    # Set tick parameters to font size 14
    ax.tick_params(axis='both', which='major', labelsize=14)

    # Set axis limits
    ax.set_xlim(0, 0.3)
    ax.set_ylim(0, 0.3)

    # Create a dictionary to store plotted types for the legend
    handles = []
    
    # Plot each unique type with a different marker and color
    unique_types = df['Type'].unique()
    for i, t in enumerate(unique_types):
        df_type = df[df['Type'] == t]
        scatter = ax.scatter(df_type['Predicted'], df_type['Observed'], label=t.replace('_', ' '), 
                             marker=markers[i % len(markers)], color=colors(i), s=100)
        handles.append(scatter)

    # Optionally, add a reference line y=x to compare
    ax.plot([0, 0.3], [0, 0.3], color='red', linestyle='--', label='Observed=Predicted')

    # Add a legend using the handles for each type, no title for the legend
    ax.legend(handles=handles, fontsize=14, loc='lower right')

    # Add grid for better readability
    ax.grid(True)

    # Adjust layout and save the figure
    plt.tight_layout()
    plt.savefig(f'{filename}.png', dpi=300)
    plt.show()

# Plot and save the scatter plot with different markers and legend
plot_scatter_with_legend(df_filtered, '', 'scatter_plot_diverged')

