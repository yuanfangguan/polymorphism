import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

import glob

all_base=glob.glob("../base_excludeN_5fold_512/auc.fold.*")
group_A_values=[]
values=[]
for the_file in all_base:
    if ('NN' in the_file):
      #  pass
   # else:
        FILE=open(the_file,'r')
        line=FILE.readline()
        val=float(line.strip())
        group_A_values.append(val)
        values.append(val)
        FILE.close()

all_base=glob.glob("../base_excludeN_5fold_512/train.auc.fold.*")
for the_file in all_base:
    if ('NN' in the_file):
   #     pass
  #  else:
        FILE=open(the_file,'r')
        line=FILE.readline()
        val=float(line.strip())
        group_A_values.append(val)
        values.append(val)
        FILE.close()


all_base=glob.glob("../base_excludeN_5fold_1024/auc.fold.*")
group_B_values=[]
for the_file in all_base:
    if ('NN' in the_file):
        FILE=open(the_file,'r')
        line=FILE.readline()
        val=float(line.strip())
        group_B_values.append(val)
        values.append(val)
        FILE.close()

all_base=glob.glob("../base_excludeN_5fold_1024/train.auc.fold.*")
for the_file in all_base:
    if ('NN' in the_file):
        FILE=open(the_file,'r')
        line=FILE.readline()
        val=float(line.strip())
        group_B_values.append(val)
        values.append(val)
        FILE.close()


all_base=glob.glob("../base_excludeN_5fold/auc.fold.*")
group_C_values=[]
for the_file in all_base:
    if ('NN' in the_file):
        FILE=open(the_file,'r')
        line=FILE.readline()
        val=float(line.strip())
        group_C_values.append(val)
        values.append(val)
        FILE.close()

all_base=glob.glob("../base_excludeN_5fold/train.auc.fold.*")
for the_file in all_base:
    if ('NN' in the_file):
        FILE=open(the_file,'r')
        line=FILE.readline()
        val=float(line.strip())
        group_C_values.append(val)
        values.append(val)
        FILE.close()


all_base=glob.glob("../base_excludeN_5fold_4096/auc.fold.*")
group_D_values=[]
for the_file in all_base:
    if ('NN' in the_file):
        FILE=open(the_file,'r')
        line=FILE.readline()
        val=float(line.strip())
        group_D_values.append(val)
        values.append(val)
        FILE.close()

all_base=glob.glob("../base_excludeN_5fold_4096/train.auc.fold.*")
for the_file in all_base:
    if ('NN' in the_file):
        FILE=open(the_file,'r')
        line=FILE.readline()
        val=float(line.strip())
        group_D_values.append(val)
        values.append(val)
        FILE.close()

data = {
    'Group': ['512 bp'] * 10 + ['1024 bp'] * 10+ ['2048 bp']*10+ ['4096 bp']*10,
    'Set': ['Test'] * 5 + ['Train'] * 5 + ['Test'] * 5 + ['Train'] * 5 +['Test']*5+['Train']*5+['Test']*5+['Train']*5,
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
plt.ylabel('AUC per nucleotide')

plt.legend(title='', loc='lower right', bbox_to_anchor=(1, 0), frameon=False)


# Set x and y limits if needed
plt.ylim(0.68, 0.73)

# Save the plot as a PNG file
plt.savefig('chunk_size_auc_pernucleotide.png', dpi=300)

# Show the plot
plt.show()

