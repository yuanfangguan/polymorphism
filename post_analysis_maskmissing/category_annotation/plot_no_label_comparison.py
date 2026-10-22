import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv('no_label_comparison.txt', sep='\t')
print(df.columns)

# Sort the DataFrame based on the Observed values
df_sorted = df.sort_values(by='Observed')

# Split the data into two halves
mid_index = len(df_sorted) // 2
df_first_half = df_sorted[:mid_index]
df_second_half = df_sorted[mid_index:]

# Function to plot the data
def plot_data(df, title, filename):
    fig, ax = plt.subplots()
    index = range(len(df))
    bar_width = 0.45
    opacity = 0.67

    rects1 = ax.bar(index, df['Observed'], bar_width, alpha=opacity, color='red', label='Observed SNP ratio')
    rects2 = ax.bar([p + bar_width for p in index], df['Annotation masked model'], bar_width, alpha=opacity, color='blue', label='Annotation-masked model')

    #ax.set_xlabel('Type', fontsize=12)
   # ax.set_ylabel('Values', fontsize=12)
    ax.set_title(title, fontsize=12)
    ax.set_xticks([p + bar_width / 2 for p in index])
    ax.set_xticklabels(df['Type'], rotation=45, ha='right', fontsize=12)
    ax.legend(fontsize=12)

    # Dynamically adjust the y-axis upper limit
    max_value = max(df['Observed'].max(), df['Annotation masked model'].max())
    plt.ylim(0, 1.25)  # Increase upper limit by 15%

    plt.tight_layout()
    plt.savefig(f'{filename}.png', dpi=300)  # Save the figure with higher DPI

# Plot and save the two halves
plot_data(df_first_half, '', 'first_half_comparison_nolabel')
plot_data(df_second_half, '', 'second_half_comparison_nolabel')

