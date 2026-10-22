import os
import numpy as np
import pandas as pd
from glob import glob
from sklearn.metrics import roc_auc_score

# ---------------------------------------
# 1) Load Label File
# ---------------------------------------

LABEL_FILE = "human.denovo.mutation.dataset1"  # name of your label file

# Example columns (customize parsing as needed):
# chr1    619713  C   T   Proband-1   ...
# chr1    786760  C   G   Proband-2   ...
# ...
#
# We'll assume columns:
#   0: chromosome
#   1: position
#
# And we label certain positions as 1, all others default to 0.
import pandas as pd


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

# We’ll collect per-chrom results here
auc_per_chrom = {}

# Optionally store all genome-wide labels/preds for an overall AUC
all_labels = []
all_preds = []

for npy_path in npy_files:
    filename = os.path.basename(npy_path)  
    # e.g. "NC_000001.11.pred.npy"
    refseq_id = filename.split(".pred.npy")[0]  
    # e.g. "NC_000001.11"
    
    if refseq_id not in map_chr:
        print(f"Warning: {refseq_id} not in map_chr mapping, skipping.")
        continue

    chrom_name = map_chr[refseq_id]  # e.g. "chr1"

    # Load prediction array
    # We'll assume index i in this array corresponds to position (i+1) on that chromosome
    # (adjust if your indexing is different)
    preds = np.load(npy_path)  # shape: (length_of_chrom, )

    # Construct the label array for this chromosome
    # By default, everything is 0
    labels = np.zeros(len(preds), dtype=int)

    # For positions that appear in label_dict, set to 1
    # (Again, watch out for off-by-one: if your positions start at 0, you may need to adjust.)
    for i in range(len(preds)):
        pos = i + 1  # if your positions are 1-based
        if (chrom_name, pos) in label_dict:
            labels[i] = 1

    # Compute AUC for this chromosome
    # Make sure you have at least 1 positive label and 1 negative label
    # or else roc_auc_score will fail.
    if np.sum(labels) == 0 or np.sum(labels) == len(labels):
        print(f"[{chrom_name}] Skipping AUC (all labels are 0 or all 1).")
        continue

    chrom_auc = roc_auc_score(labels, preds)
    auc_per_chrom[chrom_name] = chrom_auc

    # Optionally add to genome-wide arrays
    all_labels.extend(labels)
    all_preds.extend(preds)

# ---------------------------------------
# 4) Print AUC results
# ---------------------------------------
print("Per-chromosome AUC:")
for chrom, auc_val in sorted(auc_per_chrom.items()):
    print(f"{chrom}\t{auc_val:.4f}")

# Calculate genome-wide AUC if desired
if len(set(all_labels)) > 1:  # ensure we have both 0 and 1
    overall_auc = roc_auc_score(all_labels, all_preds)
    print(f"Genome-wide AUC: {overall_auc:.4f}")
else:
    print("Not enough variation in labels across the genome to compute a genome-wide AUC.")

