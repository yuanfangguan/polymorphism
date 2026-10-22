import pandas as pd
from functools import reduce
from pathlib import Path

files = [
    "annotation_auc.tsv",
    "annotation_auc_maskfunctional.tsv",
    "annotation_auc_maskfunctional_repeat.tsv",
    "annotation_auc_masklabel.tsv",
]

label_map = {
    "annotation_auc.tsv": "auc_base",
    "annotation_auc_maskfunctional.tsv": "auc_maskfunctional",
    "annotation_auc_maskfunctional_repeat.tsv": "auc_maskfunctional_repeat",
    "annotation_auc_masklabel.tsv": "auc_masklabel",
}

dfs = []
for f in files:
    fp = Path(f)
    if not fp.exists():
        raise FileNotFoundError(f"Missing file: {f}")
    df = pd.read_csv(fp, sep="\t", dtype={"annotation_type": str}, keep_default_na=True)
    # strip whitespace from column names just in case
    df.columns = df.columns.str.strip()
    # keep necessary cols if they exist
    needed = [c for c in ["annotation_type", "auc", "n_total", "n_pos", "n_neg"] if c in df.columns]
    df = df[needed].copy()
    if "auc" not in df.columns or "annotation_type" not in df.columns:
        raise ValueError(f"`annotation_type` and/or `auc` missing in {f}")
    df = df.rename(columns={"auc": label_map[f]})
    dfs.append(df)

# 1) Merge AUCs on annotation_type (outer to keep all)
merged_auc = reduce(
    lambda l, r: pd.merge(l, r.drop(columns=[c for c in ["n_total", "n_pos", "n_neg"] if c in r.columns]),
                          on="annotation_type", how="outer"),
    dfs
)

# 2) Choose a counts donor (first df that actually has all three)
counts_df = None
for df in dfs:
    if all(c in df.columns for c in ["n_total", "n_pos", "n_neg"]):
        counts_df = df[["annotation_type", "n_total", "n_pos", "n_neg"]].drop_duplicates("annotation_type")
        break

# 3) Merge counts back if available
if counts_df is not None:
    merged = merged_auc.merge(counts_df, on="annotation_type", how="left")
else:
    print("WARNING: No file contained all of n_total, n_pos, n_neg. Proceeding without counts.")
    merged = merged_auc

# 4) Nice ordering without KeyError
auc_cols = [label_map[f] for f in files if label_map[f] in merged.columns]
front = ["annotation_type"]
count_cols_present = [c for c in ["n_total", "n_pos", "n_neg"] if c in merged.columns]
col_order = front + count_cols_present + auc_cols
merged = merged[col_order]

# 5) Optional: sort by base AUC if present
sort_col = "auc_base" if "auc_base" in merged.columns else auc_cols[0]
merged = merged.sort_values(by=sort_col, ascending=False, na_position="last")

# 6) Save
out_path = "annotation_auc_comparison.tsv"
merged.to_csv(out_path, sep="\t", index=False)

print(f"Saved: {out_path}")
print("Columns:", list(merged.columns))
print(merged.head(20).to_string(index=False))

