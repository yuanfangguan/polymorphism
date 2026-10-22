import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import glob
plt.rcParams.update({'font.size': 14})

# Initialize lists for values in each group
group_A_train, group_A_test = [], []
group_B_train, group_B_test = [], []
group_C_train, group_C_test = [], []
group_D_train, group_D_test = [], []

# Load values for each group
def load_values(file_pattern, test_list, train_list):
    files = glob.glob(file_pattern)
    files = [f for f in glob.glob(file_pattern) if 'NNwise' in f]
    for idx, file in enumerate(files):
        with open(file, 'r') as f:
            line = f.readline()
            val = float(line.strip().split(' ')[1])

            if ('train' in file):
                train_list.append(val)
            else:
                test_list.append(val)

load_values("../base_excludeN_5fold_512/*cor.fold.*", group_A_test, group_A_train)
load_values("../base_excludeN_5fold_1024/*cor.fold.*", group_B_test, group_B_train)
load_values("../base_excludeN_5fold/*cor.fold.*", group_C_test, group_C_train)
load_values("../base_excludeN_5fold_4096/*cor.fold.*", group_D_test, group_D_train)

# Combine data for plotting
data = {
    'Train': group_A_train + group_B_train + group_C_train + group_D_train,
    'Test': group_A_test + group_B_test + group_C_test + group_D_test,
    'Group': ['512 bp'] * len(group_A_train) + ['1024 bp'] * len(group_B_train) +
             ['2048 bp'] * len(group_C_train) + ['4096 bp'] * len(group_D_train)
}

# Convert to DataFrame
df = pd.DataFrame(data)

# Define markers for each group
markers = {'512 bp': 'o', '1024 bp': 's', '2048 bp': 'D', '4096 bp': '^'}

# Plot scatter with different markers for each group
plt.figure(figsize=(6.5, 6))
for group, marker in markers.items():
    subset = df[df['Group'] == group]
    plt.scatter(subset['Train'], subset['Test'], label=group, marker=marker, s=50)

# Add labels and title
plt.xlabel('Train')
plt.ylabel('Test')
plt.legend(title='Group')
plt.title('Correlation per nucleotide by chunk size')
plt.ylim(0.25, 0.35)
plt.xlim(0.25, 0.35)
plt.plot([0, 1], [0, 1], color='black', linestyle='--', linewidth=1, label='y = x')

# Save and show the plot
plt.savefig('chunk_cor_nucloetide_scatter.png', dpi=300)
plt.show()

