import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

import glob

all_base=glob.glob("../base_excludeN_5fold_512/cor.fold.*")
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

all_base=glob.glob("../base_excludeN_5fold_maskfunctional/cor.fold.*")
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

all_base=glob.glob("../base_excludeN_5fold_maskfunctionalonlyintrain/cor.fold.*")
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


all_base=glob.glob("../base_excludeN_5fold_masklabel/cor.fold.*")
group_C_values=[]
for the_file in all_base:
    if ('NN' in the_file):
        pass
    else:
        FILE=open(the_file,'r')
        line=FILE.readline()
        val=float(line.strip().split(' ')[1])
        group_C_values.append(val)
        values.append(val)
        FILE.close()

all_base=glob.glob("../base_excludeN_5fold_masklabelonlyintrain/cor.fold.*")
for the_file in all_base:
    if ('NN' in the_file):
        pass
    else:
        FILE=open(the_file,'r')
        line=FILE.readline()
        val=float(line.strip().split(' ')[1])
        group_C_values.append(val)
        values.append(val)
        FILE.close()


data = {
    'Group': ['Original'] * 5 + ['Mask functional'] * 10+ ['Mask all annotated']*10,
    'Set': ['Original'] * 5 + ['Masked'] * 5 + ['Original'] * 5 + ['Masked'] * 5 +['Original']*5,
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
plt.ylabel('Correlation')

plt.legend(title='', loc='lower right', bbox_to_anchor=(1, 0), frameon=False)


# Set x and y limits if needed
plt.ylim(0.2, 0.8)

# Save the plot as a PNG file
plt.savefig('mask_cor.png', dpi=300)

# Show the plot
plt.show()

