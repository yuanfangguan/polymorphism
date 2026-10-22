from Bio import motifs
from Bio.Seq import Seq
from collections import defaultdict
import scipy.stats as stats
import random

# Step 1: Load sequences from the top_100000_windows.txt file
def load_sequences(input_file):
    sequences = []
    with open(input_file, 'r') as f:
        for line in f:
            # Extract the sequence from each line
            try:
                _, _, _, sequence = line.strip().split('\t')
                sequences.append(Seq(sequence))
            except ValueError:
                print(f"Skipping invalid line: {line.strip()}")

    print(f"Loaded {len(sequences)} sequences.")
    return sequences

# Step 2: Detect sub-motifs from the sequences
def detect_motifs(sequences, min_motif_len=3, max_motif_len=8):
    detected_motifs = set()

    # Search for motifs of different lengths
    for seq in sequences:
        for length in range(min_motif_len, max_motif_len + 1):
            for i in range(len(seq) - length + 1):
                motif = str(seq[i:i + length])
                detected_motifs.add(motif)

    print(f"Detected {len(detected_motifs)} motifs.")
    return detected_motifs

# Step 3: Count occurrences of each motif
def count_motif_occurrences(sequences, motifs):
    motif_counts = defaultdict(int)
    total_bases = sum(len(seq) for seq in sequences)

    # Count the occurrences of each motif
    for motif in motifs:
        for seq in sequences:
            motif_counts[motif] += seq.count(motif)

    return motif_counts, total_bases

# Step 4: Generate random sequences as background
def generate_random_sequences(num_sequences, seq_length):
    nucleotides = ['A', 'T', 'C', 'G']
    return [Seq(''.join(random.choices(nucleotides, k=seq_length))) for _ in range(num_sequences)]

# Step 5: Perform enrichment analysis
def motif_enrichment_analysis(sequences, motif_counts, total_bases):
    enrichment_results = []

    # Generate random background sequences
    background_sequences = generate_random_sequences(len(sequences), len(sequences[0]))

    # Count occurrences of motifs in the background
    background_counts, background_total_bases = count_motif_occurrences(background_sequences, motif_counts.keys())

    # Perform enrichment analysis using a binomial test
    for motif, observed_count in motif_counts.items():
        expected_count = (background_counts[motif] / background_total_bases) * total_bases

        if observed_count > 0:
            # Binomial test for enrichment
            p_value = stats.binom_test(observed_count, total_bases, expected_count / total_bases, alternative="greater")
            enrichment_results.append((motif, observed_count, expected_count, p_value))

    # Sort by p-value in ascending order
    enrichment_results.sort(key=lambda x: x[-1])
    return enrichment_results

# Main execution
if __name__ == "__main__":
    # Input file and parameters
    input_file = "top_100000_windows.txt"
    min_motif_len = 5
    max_motif_len = 10

    # Step 1: Load sequences from the input file
    sequences = load_sequences(input_file)

    # Step 2: Detect motifs from the sequences
    detected_motifs = detect_motifs(sequences, min_motif_len, max_motif_len)

    # Step 3: Count occurrences of each detected motif
    motif_counts, total_bases = count_motif_occurrences(sequences, detected_motifs)

    # Step 4: Perform enrichment analysis
    enrichment_results = motif_enrichment_analysis(sequences, motif_counts, total_bases)

    # Step 5: Print enrichment results
    print("\nEnrichment Analysis Results:")
    print(f"{'Motif':<10}{'Observed':<10}{'Expected':<10}{'P-value':<10}")
    for motif, observed, expected, p_value in enrichment_results:
        print(f"{motif:<10}{observed:<10}{expected:<10.2f}{p_value:<10.5f}")

