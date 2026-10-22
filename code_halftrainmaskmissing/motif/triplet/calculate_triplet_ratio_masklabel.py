import numpy as np
import os
from glob import glob
from Bio import SeqIO
from itertools import product

# -----------------------------------------------------------
# 1) Define paths and load the reference genome
# -----------------------------------------------------------
FASTA_FILE = "../../../../references/GRCh38_latest_genomic.fna"
PREDICTIONS_DIR = "../../ratio_masklabel"

# Create a dict: { "NC_000001.11": SeqRecord, "NC_000002.12": SeqRecord, ... }
# Using BioPython's SeqIO
print("Loading reference genome...")
ref_dict = SeqIO.to_dict(SeqIO.parse(FASTA_FILE, "fasta"))

# -----------------------------------------------------------
# 2) Get all possible triplets (64 combinations of A, C, G, T)
# -----------------------------------------------------------
nucleotides = ["A", "C", "G", "T"]
all_triplets = ["".join(p) for p in product(nucleotides, repeat=3)]

# We'll accumulate sums and counts
triplet_sums = {trip: 0.0 for trip in all_triplets}
triplet_counts = {trip: 0 for trip in all_triplets}

# -----------------------------------------------------------
# 3) Loop over all .npy prediction files
# -----------------------------------------------------------
npy_files = sorted(glob(os.path.join(PREDICTIONS_DIR, "*.npy")))
print(f"Found {len(npy_files)} prediction files in {PREDICTIONS_DIR}.")

for npy_path in npy_files:
    filename = os.path.basename(npy_path)
    # Example: "NC_000001.11.pred.npy" → refseq_id = "NC_000001.11"
    refseq_id = filename.replace(".pred.npy", "")
    refseq_id = refseq_id.replace(".npy", "")
    refseq_id = refseq_id.replace("processed_", "")

    # Make sure this chromosome/contig is in the reference dictionary
    if refseq_id not in ref_dict:
        print(f"Warning: {refseq_id} not found in reference FASTA. Skipping.")
        continue

    # Load the predictions into a numpy array
    # We'll assume preds[i] corresponds to the base at position i (0-based)
    preds = np.load(npy_path)
    chrom_seq = str(ref_dict[refseq_id].seq)  # The reference sequence as a string
    seq_len = len(chrom_seq)

    # Sanity check: the length of preds should match the sequence length
    if len(preds) != seq_len:
        print(f"Warning: length mismatch in {refseq_id}. "
              f"Sequence length={seq_len}, preds length={len(preds)}. Truncating to min length.")
        min_len = min(seq_len, len(preds))
    else:
        min_len = seq_len

    # -------------------------------------------------------
    # 4) Slide over the sequence in 3-mer (triplet) windows
    # -------------------------------------------------------
    # We'll iterate from 0 to min_len-3 so we can extract [i, i+1, i+2]
    # Middle position is i+1
    for i in range(min_len - 2):
        trip = chrom_seq[i : i+3].upper()  # 3 bases

        # Skip if it contains Ns
        if 'N' in trip:
            continue

        # The middle position in the triplet
        mid = i + 1

        # Get the predicted value at that middle position
        mid_val = preds[mid]

        # Accumulate in the dictionaries
        if trip in triplet_sums:  # Only if it's one of the canonical triplets
            triplet_sums[trip] += mid_val
            triplet_counts[trip] += 1

    print(f"Processed {refseq_id}: length={min_len} total triplets used.")

# -----------------------------------------------------------
# 5) Compute the average predicted value per triplet
# -----------------------------------------------------------
triplet_means = {}
for trip in all_triplets:
    count = triplet_counts[trip]
    if count > 0:
        avg_val = triplet_sums[trip] / count
        triplet_means[trip] = avg_val
    else:
        triplet_means[trip] = float('nan')  # or 0.0, or whatever you prefer

# -----------------------------------------------------------
# 6) Output or print results
# -----------------------------------------------------------
print("\nTriplet\tCount\tAveragePrediction")
for trip in sorted(all_triplets):
    print(f"{trip}\t{triplet_counts[trip]}\t{triplet_means[trip]:.6f}")

