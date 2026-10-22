import os
import numpy as np
from Bio import SeqIO
import random  # To shuffle the sequence
import sys

# Paths to genome file and prediction directories
genome_fasta = "../../../../references/GRCh38_latest_genomic.fna"
window_size = int(sys.argv[1])
step_size = int(int(sys.argv[1]) / 2)
top_k = 10000  # Number of random windows to retain per chromosome

# Extract genome sequence and apply masking, then generate random sequences
def get_random_windows(genome_fasta, map_genome):
    all_random_windows = []

    # Loop over each chromosome sequence in the FASTA file
    for record in SeqIO.parse(genome_fasta, "fasta"):
        chrom_name = record.id
        print(f"Processing {chrom_name}...")

        # Extract the sequence and ensure it's uppercase
        sequence = str(record.seq).upper()

        # Check if the chromosome has a corresponding masked genome
        if chrom_name not in map_genome:
            print(f"Masked data for {chrom_name} not found. Skipping.")
            continue

        random_windows = []

        # Iterate over the sequence with a sliding window
        for i in range(0, len(sequence) - window_size + 1, step_size):
            # Check if the current window overlaps with a masked region
            mask_window = map_genome[chrom_name][i:i + window_size, :]

            # If any base in the window is masked (all zeros), skip this window
            if np.any(np.all(mask_window == [0, 0, 0, 0], axis=1)):
                continue

            # Get the sequence for the current window
            window_seq = sequence[i:i + window_size]

            # Shuffle the sequence randomly while keeping the composition
            window_seq_list = list(window_seq)
            random.shuffle(window_seq_list)
            random_seq = ''.join(window_seq_list)

            # Store the window information
            random_windows.append((chrom_name, i, random_seq))

            # Stop if we have collected the top_k windows for this chromosome
            if len(random_windows) >= top_k:
                break

        # Add windows for this chromosome to the main list
        all_random_windows.extend(random_windows)

    return all_random_windows


# Run the process and save the random windows to a file
if __name__ == "__main__":
    GENOME = open('../../../data/GRCh38_latest_genomic.fna', 'r')
    map_ACTG = {'A': [1, 0, 0, 0], 'C': [0, 1, 0, 0], 'T': [0, 0, 1, 0], 'G': [0, 0, 0, 1]}
    map_genome = {}
    
    for line in GENOME:
        line = line.strip()
        if '>' in line:
            table = line.split(' ')
            table = table[0].split('>')
            the_chr = table[1]
            if 'NC_' in the_chr:
                map_genome[the_chr] = []
        else:
            line = line.upper()
            for aaa in line:
                if 'NC_' in the_chr:
                    try:
                        map_genome[the_chr].append(map_ACTG[aaa])
                    except KeyError:
                        map_genome[the_chr].append([0, 0, 0, 0])
    GENOME.close()

    TYPE = open('train_type.txt', 'r')
    type_map = {line.strip(): i for i, line in enumerate(TYPE)}
    TYPE.close()

    all_chr = map_genome.keys()
    for the_chr in all_chr:
        map_genome[the_chr] = np.asarray(map_genome[the_chr])

    print("start mask")
    GENOME = open('../../../../references/GRCh38_latest_genomic.gff', 'r')
    for line in GENOME:
        if '#' not in line:
            line = line.strip()
            table = line.split('\t')
            if table[1] != 'RefSeq' and table[0] in map_genome:
                try:
                    map_genome[table[0]][(int(table[3]) - 1) - 512:int(table[4]) + 512, :] = 0
                    print('masked')
                except IndexError:
                    pass
    GENOME.close()

    # Generate random windows per chromosome instead of ranking by prediction values
    random_windows = get_random_windows(genome_fasta, map_genome)

    # Save the random windows to a file
    output_file = f"random_{top_k}_{sys.argv[1]}_windows.txt"
    with open(output_file, "w") as f:
        for chrom, pos, seq in random_windows:
            f.write(f"{chrom}\t{pos}\t{seq}\n")

    print(f"Random windows saved to {output_file}")

