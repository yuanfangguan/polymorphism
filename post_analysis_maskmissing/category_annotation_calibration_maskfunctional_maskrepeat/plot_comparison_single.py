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
    rects2 = ax.bar([p + bar_width for p in index], df['Predicted'], bar_width, alpha=opacity, color='blue', label='Predicted')

    #ax.set_xlabel('Type', fontsize=14)
    #ax.set_ylabel('Values', fontsize=14)
    ax.set_title(title, fontsize=14)
    ax.set_ylabel('Frequency of mutation\n', fontsize=14)

    # Replace underscores with spaces in x-axis labels
    df['Type'] = df['Type'].str.replace('_', ' ')
    ax.set_xticks([p + bar_width / 2 for p in index])
    ax.set_xticklabels(df['Type'], rotation=45, ha='right', fontsize=14)
    ax.legend(fontsize=14, loc='upper left')
    ax.tick_params(axis='y', labelsize=14)


    # Dynamically adjust the y-axis upper limit
    max_value = max(df['Observed'].max(), df['Predicted'].max())
    plt.ylim(0, max_value * 1.1)  # Increase upper limit by 25%

    plt.tight_layout()
    plt.savefig(f'{filename}.png', dpi=300)
    plt.show()

# Plot and save the combined data
plot_combined_data(df_combined, '', 'no_label_combined_comparison')

