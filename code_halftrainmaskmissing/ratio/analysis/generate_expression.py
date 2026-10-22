import pysam
import numpy as np
import os

# Paths to your BAM file and its index
bam_file = "GCF_000001405.40-RS_2023_03_knownrefseq_alns.bam"
bam_index = "GCF_000001405.40-RS_2023_03_knownrefseq_alns.bam.bai"

# Create output directory for npy files
output_dir = "chromosome_npy_files"
os.makedirs(output_dir, exist_ok=True)

# Open the BAM file using pysam
bam = pysam.AlignmentFile(bam_file, "rb")

# Step through each chromosome in the BAM file
for chromosome in bam.references:
    chromosome_length = bam.get_reference_length(chromosome)
    
    # Initialize an array to store the coverage for each position
    coverage = np.zeros(chromosome_length, dtype=np.int32)
    
    # Fetch the pileup (coverage) for the chromosome
    for pileupcolumn in bam.pileup(chromosome):
        position = pileupcolumn.reference_pos  # 0-based position
        coverage[position] = pileupcolumn.n  # Number of reads covering this position
    
    # Save the coverage array as an npy file for the chromosome
    npy_file = os.path.join(output_dir, f"{chromosome}_coverage.npy")
    np.save(npy_file, coverage)
    print(f"Saved coverage for {chromosome} to {npy_file}")

# Close the BAM file
bam.close()

