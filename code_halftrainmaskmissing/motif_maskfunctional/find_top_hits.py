import os
import numpy as np
from Bio import SeqIO
import heapq  # For maintaining the top 100,000 windows efficiently

# Paths to genome file and prediction directories
genome_fasta = "../../../references/GRCh38_latest_genomic.fna"
prediction_dir_1 = "../../calibrated_predictions_maskfunctional/"
window_size = 128
step_size = 64
top_k = 100000  # Number of top windows to retain

# Helper function to load predictions for a chromosome
def load_prediction(chrom_name):
    pred_file = f"{prediction_dir_1}/{chrom_name}.pred.npy"
    return np.load(pred_file)

# Extract genome sequence and predictions
def get_top_windows(genome_fasta):
    # Use a min-heap to store the top 100,000 windows by average prediction
    top_windows = []

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

            # Calculate the average prediction value for this window
            avg_pred = np.mean(window_pred)

            # Store the window in the min-heap (size limited to top 100,000)
            if ('N' in window_seq):
                pass
            else:
                if len(top_windows) < top_k:
                # Add the window if we have not reached the limit
                    heapq.heappush(top_windows, (avg_pred, chrom_name, i, window_seq))
                else:
                    # If the heap is full, replace the smallest element if the new one is larger
                    heapq.heappushpop(top_windows, (avg_pred, chrom_name, i, window_seq))

    return top_windows

# Run the process and save the top 100,000 windows to a file
if __name__ == "__main__":
    top_windows = get_top_windows(genome_fasta)

    # Sort the windows by average prediction value in descending order
    top_windows = sorted(top_windows, key=lambda x: x[0], reverse=True)

    # Save the top windows to a file
    output_file = "top_100000_windows.txt"
    with open(output_file, "w") as f:
        for avg_pred, chrom, pos, seq in top_windows:
            f.write(f"{chrom}\t{pos}\t{avg_pred:.4f}\t{seq}\n")

    print(f"Top 100,000 windows saved to {output_file}")

