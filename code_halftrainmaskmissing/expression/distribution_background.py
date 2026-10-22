import numpy as np
import pickle
import pyBigWig

# Chromosome name mapping
chrom_map = {
    "NC_000001.11": "chr1", "NC_000002.12": "chr2", "NC_000003.12": "chr3", "NC_000004.12": "chr4",
    "NC_000005.10": "chr5", "NC_000006.12": "chr6", "NC_000007.14": "chr7", "NC_000008.11": "chr8",
    "NC_000009.12": "chr9", "NC_000010.11": "chr10", "NC_000011.10": "chr11", "NC_000012.12": "chr12",
    "NC_000013.11": "chr13", "NC_000014.9": "chr14", "NC_000015.10": "chr15", "NC_000016.10": "chr16",
    "NC_000017.11": "chr17", "NC_000018.10": "chr18", "NC_000019.10": "chr19", "NC_000020.11": "chr20",
    "NC_000021.9": "chr21", "NC_000022.11": "chr22", "NC_000023.11": "chrX", "NC_000024.10": "chrY"
}

# Load chromosome list
def load_chromosomes(filename):
    with open(filename, 'r') as file:
        chromosomes = [line.strip() for line in file.readlines()]
    return chromosomes

# Load annotated regions (GFF or BED file, here we assume GFF format)
def load_annotations(filename):
    annotations = {}
    with open(filename, 'r') as file:
        for line in file:
            if line.startswith("#"):  # Skip comment lines
                continue

            columns = line.strip().split('\t')
            if (columns[1]=='RefSeq'):
                pass
            else:
                chr_name, feature_type, start, end = columns[0], columns[2], int(columns[3]), int(columns[4])
                chr_name=chrom_map.get(chr_name, chr_name)
                if chr_name not in annotations:
                    annotations[chr_name] = []
                annotations[chr_name].append((start, end))
    return annotations

# Check if a region is within ±500 of any annotated region
def is_near_annotation(start, end, annotations, chromosome, distance=500):

    if chromosome not in annotations:

        print(chromosome)
        return False
    for ann_start, ann_end in annotations[chromosome]:
        if (start < ann_end + distance) and (end > ann_start - distance):
            return True
    return False

# Parameters
window_size = 64
num_samples_per_chr = 5000  # Limiting to 500 samples per chromosome
output_file = "mean_expression_10000_regions.txt"

# Load annotations
annotations = load_annotations("../../../references/GRCh38_latest_genomic.gff")  # Replace with your annotation file path

# Load chromosome list
chr_list = load_chromosomes('train.chr')

# Load bigWig file
bw = pyBigWig.open("RNAseq_default.bw")

# Processing chromosome-by-chromosome
with open(output_file, "w") as f:
    for the_chr in chr_list:
        mapped_chr = chrom_map.get(the_chr, the_chr)  # Convert to UCSC-style name
        try:
            # Load data and label for the current chromosome
            with open(f'../../data_masklabel/{the_chr}.pkl', "rb") as data_file:
                data = pickle.load(data_file)

            with open(f'../../data/{the_chr}_label.pkl', "rb") as label_file:
                label_data = pickle.load(label_file)
        except FileNotFoundError:
            print(f"Skipping {the_chr}: Data file not found.")
            continue

        # Find valid positions where label == 0
        valid_positions = np.where(label_data == 0)[0]
        if len(valid_positions) == 0:
            print(f"Skipping {mapped_chr}: No valid positions found.")
            continue

        # Randomly sample up to `num_samples_per_chr`
        sampled_positions = np.random.choice(valid_positions, size=min(num_samples_per_chr, len(valid_positions)), replace=False)

        for start_pos in sampled_positions:
            end_pos = start_pos + window_size
#            print(start_pos, end_pos)
            

            # Skip regions near annotations (±500 bp)
            if is_near_annotation(start_pos, end_pos, annotations, mapped_chr):
                continue

            try:
                # Fetch all values for the region
                values = bw.values(mapped_chr, start_pos, end_pos, numpy=True)

                # Check if any value in the region is NaN
                mean_expr = np.nanmean(values)
                f.write(f"{mean_expr}\n")  # Directly write to file

            except RuntimeError as e:
                print(f"Skipping region {mapped_chr}:{start_pos}-{end_pos} due to error.")
                print(e)

bw.close()
print(f"Saved mean expression values to {output_file}")

