import re
import glob
import os

# Define the input and output files
all_files = glob.glob('*.viennafold')

# Regular expression pattern to match floating-point numbers enclosed in parentheses
energy_pattern = re.compile(r'\((-?\d+\.\d+)\)')

for the_file in all_files:
    input_file = the_file  # ViennaRNA output filename
    output_file = the_file + '.summary'

    # List to store all extracted energy values
    all_values = []

    # Process the file line by line
    with open(input_file, 'r') as file:
        for line in file:
            # Find the energy value in parentheses at the end of the line
            match = energy_pattern.search(line)
            if match:
                # Extract the numeric part from the match and convert to float
                energy_value = float(match.group(1))  # Extract the number inside the parentheses
                all_values.append(energy_value)

    # Write the extracted values to an output file if any were found
    if all_values:
        with open(output_file, 'w') as OUT:
            for value in all_values:
                OUT.write(f"{value}\n")
    else:
        print(f"No energy values found in {input_file}.")

