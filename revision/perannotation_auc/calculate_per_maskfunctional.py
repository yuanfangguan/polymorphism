#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Compute per-annotation ROC-AUC (plus an 'unannotated' bucket) by combining:
- Per-chromosome predictions:  ../calibrated_predictions/NC_*.pred.npy
- Per-chromosome labels:       ../data/NC_*_label.pkl
- GFF genome annotations:      /local/disk4/gyuanfan/references/GRCh38_latest_genomic.gff

Key behavior:
- If predictions and labels differ slightly in length, we HEAD-ALIGN and TRUNCATE both
  to the shorter length (common in your setup, usually off by ~2048).
- GFF is parsed once. Intervals are stored 0-based half-open.
- Unannotated is taken as the complement of the union of all annotated intervals per chr.
- AUC is computed only if both classes are present; otherwise records 'NA'.
- Output: annotation_auc.tsv
"""

import os
import re
import sys
import pickle
import numpy as np
from collections import defaultdict
from sklearn.metrics import roc_auc_score

# ----------------------------
# Config — update paths if needed
# ----------------------------
PRED_DIR = "../../calibrated_predictions_maskfunctional"
LABEL_DIR = "../../data"
GFF_PATH = "/local/disk4/gyuanfan/references/GRCh38_latest_genomic.gff"
OUT_TSV = "annotation_auc_maskfunctional.tsv"

# Optional: limit to these annotation types (None = use all types seen in GFF)
# TYPES_WHITELIST = {"gene","exon","CDS"}  # example
TYPES_WHITELIST = None

# Optional: exclude features by source (your snippet excluded 'RefSeq')
EXCLUDE_SOURCES = {"RefSeq"}

UNANN = "unannotated"


# ----------------------------
# Helpers
# ----------------------------
def parse_gff_intervals(gff_path, types_whitelist=None, exclude_sources=None):
    """
    Build: ann_intervals[seqid][ann_type] = list of (start0, end1) 0-based half-open
    Also track the full set of types encountered.
    """
    ann_intervals = defaultdict(lambda: defaultdict(list))
    all_types = set()
    with open(gff_path, "r") as f:
        for line in f:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            seqid, source, ftype, start, end, score, strand, phase, attrs = parts

            if exclude_sources and source in exclude_sources:
                continue
            if types_whitelist is not None and ftype not in types_whitelist:
                continue

            try:
                start0 = int(start) - 1  # GFF is 1-based inclusive -> 0-based
                end1 = int(end)          # keep as half-open
                if end1 <= start0:
                    continue
            except Exception:
                continue

            ann_intervals[seqid][ftype].append((start0, end1))
            all_types.add(ftype)

    # Merge intervals per type per seqid to simplify downstream work
    def merge_intervals(intervals):
        if not intervals:
            return []
        intervals = sorted(intervals)
        merged = [intervals[0]]
        for s, e in intervals[1:]:
            ps, pe = merged[-1]
            if s <= pe:  # overlap/touch
                merged[-1] = (ps, max(pe, e))
            else:
                merged.append((s, e))
        return merged

    for seqid in ann_intervals:
        for ftype in list(ann_intervals[seqid].keys()):
            ann_intervals[seqid][ftype] = merge_intervals(ann_intervals[seqid][ftype])

    return ann_intervals, all_types


def iter_chr_files(pred_dir):
    """
    Yields (seqid, pred_path, label_path) for NC_* present in pred_dir with a matching label file.
    """
    for fn in os.listdir(pred_dir):
        m = re.match(r"^(NC_\d{6}\.\d+)\.pred\.npy$", fn)
        if not m:
            continue
        seqid = m.group(1)
        pred_path = os.path.join(pred_dir, fn)
        label_path = os.path.join(LABEL_DIR, f"{seqid}_label.pkl")
        if os.path.exists(label_path):
            yield seqid, pred_path, label_path


def slice_arrays_by_intervals(arr, intervals):
    """
    Given a 1D array and list of (start0,end1) non-overlapping sorted intervals,
    return a concatenated 1D array of arr over those intervals.
    Safe with bounds (intervals out of range are skipped).
    """
    if not intervals:
        return np.array([], dtype=arr.dtype)
    L = arr.shape[0]
    chunks = []
    for s, e in intervals:
        if s < 0:
            s = 0
        if e > L:
            e = L
        if s < e:
            chunks.append(arr[s:e])
    if not chunks:
        return np.array([], dtype=arr.dtype)
    return np.concatenate(chunks, axis=0)


def complement_intervals(intervals, L):
    """
    Given merged, sorted intervals within [0,L), return the complementary (unannotated) intervals.
    """
    if not intervals:
        return [(0, L)] if L > 0 else []
    comp = []
    prev = 0
    for s, e in intervals:
        if s > prev:
            comp.append((prev, s))
        prev = max(prev, e)
    if prev < L:
        comp.append((prev, L))
    return comp


# ----------------------------
# Main
# ----------------------------
def main():
    print("Parsing GFF...")
    ann_intervals, all_types = parse_gff_intervals(
        GFF_PATH,
        types_whitelist=TYPES_WHITELIST,
        exclude_sources=EXCLUDE_SOURCES,
    )
    if TYPES_WHITELIST is not None:
        all_types = set(TYPES_WHITELIST) & all_types  # tighten to what actually appears

    # Accumulators across chromosomes
    by_type_preds = defaultdict(list)
    by_type_labels = defaultdict(list)

    chr_count = 0
    for seqid, pred_path, label_path in iter_chr_files(PRED_DIR):
        chr_count += 1

        preds = np.load(pred_path)
        with open(label_path, "rb") as f:
            labels = pickle.load(f)
        labels = np.asarray(labels)

        # --- HEAD-ALIGN & TRUNCATE if lengths differ ---
        if labels.shape[0] != preds.shape[0]:
            L = min(preds.shape[0], labels.shape[0])
            print(f"[INFO] Length mismatch for {seqid}: preds {preds.shape[0]} vs labels {labels.shape[0]} -> trunc to {L} (head-aligned).")
            preds = preds[:L]
            labels = labels[:L]
        else:
            L = preds.shape[0]
        if L == 0:
            print(f"[WARN] Empty after truncation for {seqid}, skipping.")
            continue
        # -----------------------------------------------

        # If this chromosome has no GFF entries, treat all positions as unannotated
        seq_ann = ann_intervals.get(seqid, {})
        types_here = set(seq_ann.keys()) if seq_ann else set()
        types_to_use = (all_types & types_here) if all_types else types_here

        # Track all annotated intervals (across all types) to later derive unannotated
        all_ann_for_chr = []

        # Collect per-type
        for ftype in sorted(types_to_use):
            ivals = seq_ann[ftype]
            all_ann_for_chr.extend(ivals)
            p_sub = slice_arrays_by_intervals(preds, ivals)
            y_sub = slice_arrays_by_intervals(labels, ivals)
            if p_sub.size > 0 and y_sub.size > 0:
                by_type_preds[ftype].append(p_sub)
                by_type_labels[ftype].append(y_sub)

        # Build unannotated as complement of merged(all annotated)
        if all_ann_for_chr:
            all_ann_for_chr = sorted(all_ann_for_chr)
            merged = [all_ann_for_chr[0]]
            for s, e in all_ann_for_chr[1:]:
                ps, pe = merged[-1]
                if s <= pe:
                    merged[-1] = (ps, max(pe, e))
                else:
                    merged.append((s, e))
            comp = complement_intervals(merged, L)
        else:
            comp = [(0, L)]

        # Slice unannotated
        p_un = slice_arrays_by_intervals(preds, comp)
        y_un = slice_arrays_by_intervals(labels, comp)
        if p_un.size > 0 and y_un.size > 0:
            by_type_preds[UNANN].append(p_un)
            by_type_labels[UNANN].append(y_un)

    if chr_count == 0:
        print("No chromosomes found. Check PRED_DIR and LABEL_DIR.")
        sys.exit(1)

    # Compute AUC per type
    results = []
    all_keys = sorted(by_type_preds.keys())
    for ftype in all_keys:
        p = np.concatenate(by_type_preds[ftype]) if by_type_preds[ftype] else np.array([])
        y = np.concatenate(by_type_labels[ftype]) if by_type_labels[ftype] else np.array([])

        n = int(y.size)
        n_pos = int(y.sum()) if n > 0 else 0
        n_neg = int(n - n_pos)

        if n_pos > 0 and n_neg > 0:
            try:
                auc = roc_auc_score(y, p)
            except Exception:
                auc = float("nan")
        else:
            auc = float("nan")

        results.append((ftype, auc, n, n_pos, n_neg))

    # Write TSV
    with open(OUT_TSV, "w") as f:
        f.write("annotation_type\tauc\tn_total\tn_pos\tn_neg\n")
        # Put unannotated first, then alphabetical
        for ftype, auc, n, n_pos, n_neg in sorted(results, key=lambda x: (x[0] != UNANN, x[0])):
            auc_str = f"{auc}" if not (isinstance(auc, float) and np.isnan(auc)) else "NA"
            f.write(f"{ftype}\t{auc_str}\t{n}\t{n_pos}\t{n_neg}\n")

    print(f"Done. Wrote {OUT_TSV}")
    # Also print a brief summary
    preview = sorted(results, key=lambda x: (x[0] != UNANN, x[0]))[:20]
    for ftype, auc, n, n_pos, n_neg in preview:
        auc_s = f"{auc:.4f}" if not (isinstance(auc, float) and np.isnan(auc)) else "NA"
        print(f"{ftype:20s}  AUC={auc_s}  N={n}  pos={n_pos}  neg={n_neg}")


if __name__ == "__main__":
    main()

