import os
import numpy as np
import pandas as pd
from glob import glob
from sklearn.metrics import roc_auc_score

# ---------------------------------------
# 1) Load Label File
# ---------------------------------------

LABEL_FILE = "human.denovo.mutation.dataset1"  # name of your label file

labels_df = pd.read_csv(
    LABEL_FILE,
    sep="\t",
    header=None,    # Assuming no header in your file
    usecols=[0, 1]  # Only read columns 0 and 1
)

# Rename them so they make sense in your code
labels_df.columns = ["chrom", "pos"]

# Now create a dictionary for labeling
label_dict = {}
for _, row in labels_df.iterrows():
    c = row["chrom"]
    p = int(row["pos"])
    label_dict[(c, p)] = 1

# ---------------------------------------
# 2) Map RefSeq Names to Chromosome Names
# ---------------------------------------
map_chr = {
    "NC_000001.11": "chr1",
    "NC_000002.12": "chr2",
    "NC_000003.12": "chr3",
    "NC_000004.12": "chr4",
    "NC_000005.10": "chr5",
    "NC_000006.12": "chr6",
    "NC_000007.14": "chr7",
    "NC_000008.11": "chr8",
    "NC_000009.12": "chr9",
    "NC_000010.11": "chr10",
    "NC_000011.10": "chr11",
    "NC_000012.12": "chr12",
    "NC_000013.11": "chr13",
    "NC_000014.9":  "chr14",
    "NC_000015.10": "chr15",
    "NC_000016.10": "chr16",
    "NC_000017.11": "chr17",
    "NC_000018.10": "chr18",
    "NC_000019.10": "chr19",
    "NC_000020.11": "chr20",
    "NC_000021.9":  "chr21",
    "NC_000022.11": "chr22",
    "NC_000023.11": "chrX",
    "NC_000024.10": "chrY",
}

# ---------------------------------------
# 3) Load Predictions and Calculate AUC
# ---------------------------------------
PREDICTION_DIR = "../calibrated_predictions"
npy_files = glob(os.path.join(PREDICTION_DIR, "*.npy"))  # e.g., NC_000001.11.pred.npy, etc.

auc_per_chrom = {}

# Optionally store all genome-wide labels/preds for an overall AUC
all_labels = []
all_preds = []

# Probability to keep a negative site
NEG_KEEP_PROB = 1 / 100000  # 1 in 1,000,000

for npy_path in npy_files:
    filename = os.path.basename(npy_path)
    refseq_id = filename.replace(".pred.npy", "")
    
    if refseq_id not in map_chr:
        print(f"Warning: {refseq_id} not in map_chr mapping, skipping.")
        continue

    chrom_name = map_chr[refseq_id]  # e.g. "chr1"

    # Load prediction array
    preds = np.load(npy_path)  # shape: (length_of_chrom, )

    # Construct the label array for this chromosome
    labels = np.zeros(len(preds), dtype=int)

    # Mark the positives as 1
    # If your positions are 1-based in label_dict, then pos = i + 1
    for i in range(len(preds)):
        pos = i + 1  # if your positions are 1-based
        if (chrom_name, pos) in label_dict:
            labels[i] = 1

    # Now we subsample the negatives (label == 0)
    # Keep all positives
    mask_positives = (labels == 1)
    # For the negatives, keep a random ~1 in 1,000,000
    random_vals = np.random.rand(len(labels))
    mask_negatives_subsample = (labels == 0) & (random_vals < NEG_KEEP_PROB)

    # Combine them
    final_mask = mask_positives | mask_negatives_subsample

    subsampled_labels = labels[final_mask]
    subsampled_preds = preds[final_mask]

    # Check we have at least 1 positive and 1 negative
    if len(np.unique(subsampled_labels)) < 2:
        print(f"[{chrom_name}] Skipping AUC (subsampled data has no positive or no negative).")
        continue

    # Compute AUC on the subset
    chrom_auc = roc_auc_score(subsampled_labels, subsampled_preds)
    auc_per_chrom[chrom_name] = chrom_auc

    # Optionally add to genome-wide arrays
    all_labels.extend(subsampled_labels)
    all_preds.extend(subsampled_preds)

# ---------------------------------------
# 4) Print AUC results
# ---------------------------------------
print("Per-chromosome AUC (with negative subsampling):")
for chrom, auc_val in sorted(auc_per_chrom.items()):
    print(f"{chrom}\t{auc_val:.4f}")

# Calculate genome-wide AUC if desired
if len(set(all_labels)) > 1:  # ensure we have at least one 0 and one 1
    overall_auc = roc_auc_score(all_labels, all_preds)
    print(f"Genome-wide AUC (with negative subsampling): {overall_auc:.4f}")
else:
    print("Not enough variation in the subsampled labels to compute a genome-wide AUC.")

