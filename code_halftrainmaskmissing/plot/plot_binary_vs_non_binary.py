import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

import glob

all_base=glob.glob("../base_excludeN_5fold/auc.fold.NNwise*")
group_A_values=[]
for the_file in all_base:
    FILE=open(the_file,'r')
    line=FILE.readline()
    val=float(line.strip())
    group_A_values.append(val)
    FILE.close()

all_base=glob.glob("../base_excludeN_nonbinary_5fold/auc.fold.NNwise*")
group_B_values=[]
for the_file in all_base:
    FILE=open(the_file,'r')
    line=FILE.readline()
    val=float(line.strip())
    group_B_values.append(val)
    FILE.close()






data = {
    'Group': ['binary'] * 5 + ['continuous'] * 5,
    'Values': group_A_values + group_B_values
}

# Convert data to a DataFrame
df = pd.DataFrame(data)

# Create the box plot
plt.figure(figsize=(5, 3.5))
sns.boxplot(x='Group', y='Values', data=df)

# Add title and labels
#plt.title('Comparison of binary vs. continuous label')
plt.xlabel('')
plt.ylabel('Per-nucleotide AUC')
plt.ylim(0.6,0.7) 

# Adding text annotations to label the groups
plt.text(0, max(group_A_values) + 1, 'binary label', horizontalalignment='center', size='medium', color='black', weight='semibold')
plt.text(1, max(group_B_values) + 1, 'continuous label', horizontalalignment='center', size='medium', color='black', weight='semibold')

# Save the plot as a PNG file
plt.savefig('binary_vs_continuous_per_nucleotide.png', dpi=300)

# Show the plot
plt.show()

