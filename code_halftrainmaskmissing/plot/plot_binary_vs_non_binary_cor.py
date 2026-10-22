import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

import glob

all_base=glob.glob("../base_excludeN_5fold/cor.fold.*")
group_A_values=[]
values=[]
for the_file in all_base:
    if ('NN' in the_file):
        pass
    else:
        FILE=open(the_file,'r')
        line=FILE.readline()
        val=float(line.strip().split(' ')[1])
        group_A_values.append(val)
        values.append(val)
        FILE.close()

all_base=glob.glob("../base_excludeN_5fold/train.cor.fold.*")
for the_file in all_base:
    if ('NN' in the_file):
        pass
    else:
        FILE=open(the_file,'r')
        line=FILE.readline()
        val=float(line.strip().split(' ')[1])
        group_A_values.append(val)
        values.append(val)
        FILE.close()


all_base=glob.glob("../base_excludeN_nonbinary_5fold/cor.fold.*")
group_B_values=[]
for the_file in all_base:
    if ('NN' in the_file):
        pass
    else:
        FILE=open(the_file,'r')
        line=FILE.readline()
        val=float(line.strip().split(' ')[1])
        group_B_values.append(val)
        values.append(val)
        FILE.close()

all_base=glob.glob("../base_excludeN_nonbinary_5fold/train.cor.fold.*")
for the_file in all_base:
    if ('NN' in the_file):
        pass
    else:
        FILE=open(the_file,'r')
        line=FILE.readline()
        val=float(line.strip().split(' ')[1])
        group_B_values.append(val)
        values.append(val)
        FILE.close()

data = {
    'Group': ['binary'] * 10 + ['continuous'] * 10,
    'Set': ['Test'] * 5 + ['Train'] * 5 + ['Test'] * 5 + ['Train'] * 5,
    'Values': values
}

# Convert data to a DataFrame
df = pd.DataFrame(data)

# Create the box plot
plt.figure(figsize=(5, 3))
sns.boxplot(x='Group', y='Values', hue='Set', data=df,palette=['blue', 'red'],dodge=True)


# Add title and labels
#plt.title('Box Plot of Two Groups with Two Sets of Comparisons')
plt.xlabel('')
plt.ylabel('Correlation of 2048 chunck')
plt.legend(title='')


# Set x and y limits if needed
plt.ylim(0.2, 0.8)

# Save the plot as a PNG file
plt.savefig('binary_non_binary_cor.png', dpi=300)

# Show the plot
plt.show()

