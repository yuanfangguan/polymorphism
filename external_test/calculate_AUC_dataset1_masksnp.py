import os
import numpy as np
import pandas as pd
from glob import glob
from sklearn.metrics import roc_auc_score

# ---------------------------------------
# 1) Load Label File (de novo mutations)
# ---------------------------------------

LABEL_FILE = "human.denovo.mutation.dataset1"  # name of your label file

labels_df = pd.read_csv(
    LABEL_FILE,
    sep="\t",
    header=None,    # Assuming no header in your file
    usecols=[0, 1]  # Only read columns 0 and 1
)
labels_df.columns = ["chrom", "pos"]

# Create a dictionary of positives (de novo mutations)
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
    "NC_012920.1":  "chrM",  # if you also have the mitochondrial ref
}

# ---------------------------------------
# 3) Load Predictions, Mask SNPs, and Calculate AUC
# ---------------------------------------
PREDICTION_DIR = "../calibrated_predictions"
npy_files = glob(os.path.join(PREDICTION_DIR, "*.npy"))  # e.g., NC_000001.11.pred.npy, etc.

SNP_LABEL_DIR = "/local/disk4/gyuanfan/rate/data"  
# e.g. NC_000001.11_label.pkl, etc.

auc_per_chrom = {}
all_labels = []
all_preds = []

NEG_KEEP_PROB = 1 / 1_000_000  # 1 in 1,000,000

for npy_path in npy_files:
    filename = os.path.basename(npy_path)
    refseq_id = filename.replace(".pred.npy", "")
    
    if refseq_id not in map_chr:
        print(f"Warning: {refseq_id} not in map_chr mapping, skipping.")
        continue

    chrom_name = map_chr[refseq_id]  # e.g. "chr1"

    # -----------------------------
    # Load prediction array
    # -----------------------------
    preds = np.load(npy_path)  # shape: (length_of_chrom, )

    # -----------------------------
    # Load SNP label (binary mask)
    # This file should have the same length as preds
    # where 1 means "known SNP" and 0 means "not a known SNP."
    # -----------------------------
    snp_label_path = os.path.join(SNP_LABEL_DIR, f"{refseq_id}_label.pkl")
    if not os.path.exists(snp_label_path):
        print(f"Warning: no SNP label file found for {refseq_id}: {snp_label_path}")
        continue

    snp_label = pd.read_pickle(snp_label_path)# shape: (length_of_chrom,)
    if len(snp_label) != len(preds):
        print(f"Warning: length mismatch in SNP label for {refseq_id}")
        #continue
    
    # -----------------------------
    # Build label array for de novo sites
    # -----------------------------
    labels = np.zeros(len(preds), dtype=int)
    for i in range(len(preds)):
        pos = i + 1  # 1-based
        if (chrom_name, pos) in label_dict:
            labels[i] = 1

    # -----------------------------
    # Subsample negatives
    # Keep all positives
    # Keep 1 in 1,000,000 of the negatives
    # -----------------------------
    mask_positives = (labels == 1)
    random_vals = np.random.rand(len(labels))
    mask_negatives_subsample = (labels == 0) & (random_vals < NEG_KEEP_PROB)
    base_mask = mask_positives | mask_negatives_subsample

    # -----------------------------
    # Mask out known SNPs
    # snp_label == 1 → known SNP → exclude
    # so we keep only positions where snp_label == 0
    # -----------------------------
    mask_not_snp = (snp_label == 0)
    L_preds = len(preds)
    L_snp   = len(snp_label)
    L_min   = min(L_preds, L_snp)

    # Truncate both to the shorter length
    preds = preds[:L_min]
    snp_label = snp_label[:L_min]
    labels = labels[:L_min]  # if you also have a 'labels' array

    # Now preds, snp_label, and labels all have length = L_min
    # You can safely do:
    mask_not_snp = (snp_label == 0)
    final_mask = base_mask & mask_not_snp


    # Combine both masks

    # -----------------------------
    # Apply final mask
    # -----------------------------
    subsampled_labels = labels[final_mask]
    subsampled_preds = preds[final_mask]

    # -----------------------------
    # Compute AUC
    # -----------------------------
    unique_vals = np.unique(subsampled_labels)
    if len(unique_vals) < 2:
        print(f"[{chrom_name}] Skipping AUC (subsampled data has no positive or no negative).")
        continue

    chrom_auc = roc_auc_score(subsampled_labels, subsampled_preds)
    auc_per_chrom[chrom_name] = chrom_auc

    # Optionally add to genome-wide arrays
    all_labels.extend(subsampled_labels)
    all_preds.extend(subsampled_preds)

# ---------------------------------------
# 4) Print AUC results
# ---------------------------------------
print("Per-chromosome AUC (with negative subsampling and SNP masking):")
for chrom, auc_val in sorted(auc_per_chrom.items()):
    print(f"{chrom}\t{auc_val:.4f}")

# Calculate genome-wide AUC if desired
unique_global = np.unique(all_labels)
if len(unique_global) > 1:
    overall_auc = roc_auc_score(all_labels, all_preds)
    print(f"Genome-wide AUC (with negative subsampling & SNP masking): {overall_auc:.4f}")
else:
    print("Not enough variation in labels (after masking) to compute a genome-wide AUC.")

