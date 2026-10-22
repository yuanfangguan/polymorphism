import os
import re
import matplotlib.pyplot as plt
import seaborn as sns

# Directory containing the files
directory = './'
plt.rcParams.update({'font.size': 14})

# Define the sizes and the new categories to compare
sizes = [64, 128, 256, 512, 1024, 2048]
categories = ['lowest', 'phylop']
names = ['Top conserved regions (lowest)', 'phyloP100 (highest)']

def load_data(filename):
    """
    Parses raw .viennafold files to extract free energy values.
    Looks for the pattern ( -XX.XX) at the end of the structure line.
    """
    file_path = os.path.join(directory, filename)
    energies = []
    
    if os.path.isfile(file_path):
        with open(file_path, 'r') as f:
            for line in f:
                # Search for the energy value in parentheses at the end of dot-bracket lines
                # e.g., "....)))))) ( -7.60)" or "....)))))) (-260.50)"
                match = re.search(r'\(\s*(-?\d+\.\d+)\s*\)\s*$', line)
                if match:
                    energies.append(float(match.group(1)))
        return energies
    else:
        print(f"File {filename} not found.")
        return None

# Loop through each size and create a density plot
for size in sizes:
    plt.figure(figsize=(8, 6))
    
    for i, category in enumerate(categories):
        # Construct the filename pattern based on the category
        if category == 'lowest':
            # Assuming your original lowest files were named this way
            filename = f'lowest_100_{size}_windows.txt.viennafold'
        elif category == 'phylop':
            # Using the new phyloP files you provided
            filename = f'highest_100_{size}_phylop_windows.txt.viennafold'
            
        # If you still have .summary files, you can append '.summary' to the 
        # filenames above and revert load_data back to np.loadtxt()
        
        # Load the data
        data = load_data(filename)
        if data and len(data) > 0:
            # Plot the density curve
            sns.kdeplot(data, label=names[i], fill=True, alpha=0.3)
            
    # Set the title and labels
    plt.title(f'Segment size {size}')
    plt.xlabel('Free energy (kcal/mol)')
    plt.ylabel('Density')
    
    # Add legend
    plt.legend(loc='upper left')
    
    # Save the plot
    plt.savefig(f'density_plot_lowest_vs_phylop_{size}.png', dpi=300, bbox_inches='tight')
    plt.close() # Use close instead of show() inside loops to prevent memory leaks
