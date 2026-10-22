import matplotlib.pyplot as plt
import pandas as pd

# Assuming 'df' is your DataFrame loaded with the dataset
df = pd.read_csv('combined_comparison.txt', sep='\t')
plt.figure(figsize=(8, 6))
plt.scatter(df['Excluding annotations'], df['Whole genome model'], color='blue', alpha=0.5)
plt.title('Scatter Plot of Excluding Annotations vs. Whole Genome Model')
plt.xlabel('Excluding Annotations')
plt.ylabel('Whole Genome Model')
plt.grid(True)
plt.savefig('scatter_1.png')

