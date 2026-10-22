import os
import re
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# Directory containing the files
directory = './'
plt.rcParams.update({'font.size': 14})

# Define the sizes and the new categories to compare
sizes = [64, 128, 256, 512, 1024, 2048]
categories = ['lowest', 'phylop']
names = ['Top conserved regions (lowest)', 'phyloP100 (highest)']

# 定义两组数据的颜色：蓝色和橙色
colors = ['#1f77b4', '#ff7f0e'] 

def load_data(filename):
    """
    Parses raw .viennafold files to extract free energy values.
    """
    file_path = os.path.join(directory, filename)
    energies = []
    
    if os.path.isfile(file_path):
        with open(file_path, 'r') as f:
            for line in f:
                match = re.search(r'\(\s*(-?\d+\.\d+)\s*\)\s*$', line)
                if match:
                    energies.append(float(match.group(1)))
        return energies
    else:
        print(f"File {filename} not found.")
        return None

# Loop through each size and create a histogram plot
for size in sizes:
    plt.figure(figsize=(8, 6))
    
    loaded_datasets = []
    all_data_for_size = [] 
    
    for i, category in enumerate(categories):
        if category == 'lowest':
            filename = f'lowest_100_{size}_windows.txt.viennafold'
        elif category == 'phylop':
            filename = f'highest_100_{size}_phylop_windows.txt.viennafold'
            
        data = load_data(filename)
        if data and len(data) > 0:
            # 这里把索引 i 也存下来，方便后面取对应的颜色
            loaded_datasets.append((data, names[i], colors[i]))
            all_data_for_size.extend(data) 
            
    if not all_data_for_size:
        plt.close()
        continue
        
    min_val = min(all_data_for_size)
    max_val = max(all_data_for_size)
    common_bins = np.linspace(min_val, max_val, 30)
    
    # 3. 使用统一的 common_bins 和分配好的颜色画图
    for data, name, color in loaded_datasets:
        sns.histplot(data, label=name, bins=common_bins, color=color, alpha=0.5, stat='count', kde=False)
            
    # Set the title and labels
    plt.title(f'Segment size {size}')
    plt.xlabel('Free energy (kcal/mol)')
    plt.ylabel('Count')
    
    # Add legend
    plt.legend(loc='upper left')
    
    # Save the plot
    plt.savefig(f'histogram_lowest_vs_phylop_{size}.png', dpi=300, bbox_inches='tight')
    plt.close()
