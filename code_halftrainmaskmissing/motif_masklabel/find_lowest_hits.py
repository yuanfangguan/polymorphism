import os
import numpy as np
from Bio import SeqIO
import heapq  # For maintaining the lowest 100,000 windows efficiently

# Paths to genome file and prediction directories
genome_fasta = "../../../references/GRCh38_latest_genomic.fna"
prediction_dir_1 = "../../calibrated_predictions_masklabel/"
window_size = 128
step_size = 64
top_k = 100000  # Number of lowest windows to retain

# Helper function to load predictions for a chromosome
def load_prediction(chrom_name):
    pred_file = f"{prediction_dir_1}/{chrom_name}.pred.npy"
    return np.load(pred_file)

# Extract genome sequence and predictions
def get_lowest_windows(genome_fasta):
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

        # Iterate over the sequence with a sliding window
        for i in range(0, len(sequence) - window_size + 1, step_size):
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

            # Store the window in the max-heap (size limited to lowest 100,000)
            if ('N' in window_seq):
                pass
            else:
                if len(lowest_windows) < top_k:
                # Add the window if we have not reached the limit
                    heapq.heappush(lowest_windows, (-avg_pred, chrom_name, i, window_seq))
                else:
                # If the heap is full, replace the largest element if the new one is smaller
                    heapq.heappushpop(lowest_windows, (-avg_pred, chrom_name, i, window_seq))

    return lowest_windows

# Run the process and save the lowest 100,000 windows to a file
if __name__ == "__main__":
    lowest_windows = get_lowest_windows(genome_fasta)

    # Sort the windows by average prediction value in ascending order
    lowest_windows = sorted(lowest_windows, key=lambda x: x[0])

    # Save the lowest windows to a file
    output_file = "lowest_100000_windows.txt"
    with open(output_file, "w") as f:
        for avg_pred, chrom, pos, seq in lowest_windows:
            f.write(f"{chrom}\t{pos}\t{-avg_pred:.4f}\t{seq}\n")  # Note the negation of avg_pred to revert it

    print(f"Lowest 100,000 windows saved to {output_file}")

