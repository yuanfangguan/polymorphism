import pandas as pd
import matplotlib.pyplot as plt

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

# Plot the combined data
def plot_combined_data(df, title, filename):
    fig, ax = plt.subplots(figsize=(20, 7))
    index = range(len(df))
    bar_width = 0.35
    opacity = 0.8

    rects1 = ax.bar(index, df['Observed'], bar_width, alpha=opacity, color='red', label='Observed')
    rects2 = ax.bar([p + bar_width for p in index], df['Whole genome model'], bar_width, alpha=opacity, color='blue', label='Whole genome model')

    #ax.set_xlabel('Type', fontsize=12)
    #ax.set_ylabel('Values', fontsize=12)
    ax.set_title(title, fontsize=12)
    ax.set_xticks([p + bar_width / 2 for p in index])
    ax.set_xticklabels(df['Type'], rotation=45, ha='right', fontsize=12)
    ax.legend(fontsize=12)

    # Dynamically adjust the y-axis upper limit
    max_value = max(df['Observed'].max(), df['Whole genome model'].max())
    plt.ylim(0, max_value * 1.25)  # Increase upper limit by 25%

    plt.tight_layout()
    plt.savefig(f'{filename}.png', dpi=300)
    plt.show()

# Plot and save the combined data
plot_combined_data(df_combined, '', 'combined_comparison')

