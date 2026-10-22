import os
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# Directory containing the files
directory = './'
plt.rcParams.update({'font.size': 14})
# Define the sizes we are interested in and their corresponding categories
sizes = [64,128, 256, 512, 1024, 2048]
categories = ['lowest', 'simulated', 'random']
names= ['top conserve regions', 'simulated with same composition', 'random unannonated regions']

# Function to load data from the files
def load_data(filename):
    file_path = os.path.join(directory, filename)
    if os.path.isfile(file_path):
        return np.loadtxt(file_path)
    else:
        print(f"File {filename} not found.")
        return None

# Loop through each size and create a density plot
for size in sizes:
    plt.figure(figsize=(8, 6))
    
    i=0
    for category in categories:
        # Construct the filename pattern
        if category == 'simulated':
            filename = f'100times_simulated_{size}.txt.viennafold.summary'
        elif category == 'random':
            filename = f'{category}_10000_{size}_windows.txt.viennafold.summary'
        else:
            filename = f'{category}_100_{size}_windows.txt.viennafold.summary'

        
        # Load the data
        data = load_data(filename)
        if data is not None:
            # Plot the density curve
            sns.kdeplot(data, label=f'{names[i]}', fill=True, alpha=0.3)
        i=i+1
    
    # Set the title and labels
    plt.title(f'Segment size {size}')
    plt.xlabel('Free energy')
    plt.ylabel('Density')
    
    # Add legend
    plt.legend(loc='upper left')
    
    # Save the plot
    plt.savefig(f'density_plot_{size}.png', dpi=300)
    plt.show()

