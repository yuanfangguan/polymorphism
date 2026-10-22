import os
import numpy as np
from Bio import SeqIO
import heapq  # For maintaining the lowest 100,000 windows efficiently

# Paths to genome file and prediction directories
genome_fasta = "../../../../references/GRCh38_latest_genomic.fna"
prediction_dir_1 = "../"
window_size = 512
step_size = 256
top_k = 1000# Number of lowest windows to retain

# Helper function to load predictions for a chromosome
def load_prediction(chrom_name):
    pred_file = f"{prediction_dir_1}/processed_{chrom_name}.npy"
    return np.load(pred_file)

# Extract genome sequence and predictions
def get_lowest_windows(genome_fasta, map_genome):
    # Use a max-heap to store the lowest 100,000 windows by average prediction
    lowest_windows = []

    # Loop over each chromosome sequence in the FASTA file
    for record in SeqIO.parse(genome_fasta, "fasta"):
        chrom_name = record.id
        print(f"Processing {chrom_name}...")

        # Load the prediction values for this chromosome
        try:
            predictions = load_prediction(chrom_name)
        except FileNotFoundError:
            print(f"Prediction file for {chrom_name} not found. Skipping.")
            continue

        # Extract the sequence and ensure it's uppercase
        sequence = str(record.seq).upper()

        # Check if the chromosome has a corresponding masked genome
        if chrom_name not in map_genome:
            print(f"Masked data for {chrom_name} not found. Skipping.")
            continue

        # Iterate over the sequence with a sliding window
        for i in range(0, len(sequence) - window_size + 1, step_size):
            # Check if the current window overlaps with a masked region
            mask_window = map_genome[chrom_name][i:i + window_size, :]

            # If any base in the window is masked (all zeros), skip this window
            if np.any(np.all(mask_window == [0, 0, 0, 0], axis=1)):
                continue

            # Get the sequence and prediction values for the current window
            window_seq = sequence[i:i + window_size]
            window_pred = predictions[i:i + window_size]

            # Calculate the sum of the prediction values for this window
            window_sum = np.sum(window_pred)

            # Exclude windows where the sum is less than or equal to 0
            if window_sum <= 0:
                continue

            # Calculate the average prediction value for this window
            avg_pred = np.mean(window_pred)
            print(avg_pred)

            # Store the window in the max-heap (size limited to lowest 100,000)
            if len(lowest_windows) < top_k:
                # Add the window if we have not reached the limit
                heapq.heappush(lowest_windows, (-avg_pred, chrom_name, i, window_seq))
            else:
                # If the heap is full, replace the largest element if the new one is smaller
                heapq.heappushpop(lowest_windows, (-avg_pred, chrom_name, i, window_seq))

    return lowest_windows


# Run the process and save the lowest 100,000 windows to a file
if __name__ == "__main__":
    GENOME=open('../../../data/GRCh38_latest_genomic.fna','r')
    map_ACTG={}
    map_ACTG['A']=[1,0,0,0]
    map_ACTG['C']=[0,1,0,0]
    map_ACTG['T']=[0,0,1,0]
    map_ACTG['G']=[0,0,0,1]
    map_genome={}
    for line in GENOME:
        line=line.strip()
        if ('>' in line):
              table=line.split(' ')
              table=table[0].split('>')
              the_chr=table[1]
              if ('NC_' in the_chr):
                  map_genome[the_chr]=[]
        else:
              line=line.upper()
              for aaa in line:
                  if ('NC_' in the_chr):
                      try:
                          map_genome[the_chr].append(map_ACTG[aaa])
                      except:
                          map_genome[the_chr].append([0,0,0,0])
    GENOME.close()

    TYPE=open('train_type.txt','r')
    i=0
    type_map={}
    for line in TYPE:
        line=line.strip()
        type_map[line]=i
        i=i+1
    TYPE.close()

    all_chr=map_genome.keys()
    for the_chr in all_chr:
            map_genome[the_chr]=np.asarray(map_genome[the_chr])

    print("start mask")
    GENOME=open('../../../../references/GRCh38_latest_genomic.gff','r')
    for line in GENOME:
        if ('#' in line):
            pass
        else:
            line=line.strip()
            table=line.split('\t')
            if (table[1] == 'RefSeq'):
                pass
            else:
                print(table[0])
                if (table[0] in map_genome):
                    print(map_genome[table[0]][0:5,:])
                    try:
                        map_genome[table[0]][(int(table[3])-1)-512:int(table[4])+512,:]=0
                        print('masked')
                    except:
                        pass
    GENOME.close()

    lowest_windows = get_lowest_windows(genome_fasta, map_genome)


    # Sort the windows by average prediction value in ascending order
    lowest_windows = sorted(lowest_windows, key=lambda x: x[0])


    # Save the lowest windows to a file
    output_file = "lowest_100_512_windows.txt"
    with open(output_file, "w") as f:
        for avg_pred, chrom, pos, seq in lowest_windows:
            f.write(f"{chrom}\t{pos}\t{-avg_pred:.4f}\t{seq}\n")  # Note the negation of avg_pred to revert it

    print(f"Lowest 100,000 windows saved to {output_file}")

