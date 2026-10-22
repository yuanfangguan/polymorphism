#!/usr/bin/env python3
"""
Sanity check for ViennaRNA MFE results on window files.

Input file format (tab separated), as written by your scripts:
    chrom <TAB> pos <TAB> score <TAB> sequence

Usage:
    python sanity_check_fold.py --size 64 lowest_100_64_windows.txt highest_100_64_phylop_windows.txt
    python sanity_check_fold.py --size 128 fileA.txt fileB.txt --n 200

Optional: cross-check against the RNAfold command line binary with --cli
Optional: compare to your own energies with --energies my_energies.txt (one value per line,
          same order as the first input file)

Requires: pip install ViennaRNA numpy
"""
import argparse
import random
import subprocess
import sys

import numpy as np

try:
    import RNA
except ImportError:
    sys.exit("ViennaRNA python bindings not found. Install with: pip install ViennaRNA")


def read_windows(path, n=None):
    rows = []
    with open(path) as f:
        for ln, line in enumerate(f, 1):
            line = line.rstrip("\n")
            if not line.strip():
                continue
            parts = line.split("\t")
            if len(parts) < 4:
                print(f"  [warn] {path}:{ln} has {len(parts)} columns, expected 4")
                continue
            rows.append((parts[0], parts[1], parts[2], parts[3].strip()))
            if n and len(rows) >= n:
                break
    return rows


def check_seqs(rows, expected_len, label):
    """Check length and alphabet of every sequence."""
    bad_len, bad_chars, lower = 0, 0, 0
    lens = []
    for chrom, pos, score, seq in rows:
        lens.append(len(seq))
        if len(seq) != expected_len:
            bad_len += 1
        if any(c.islower() for c in seq):
            lower += 1
        if set(seq.upper()) - set("ACGTU"):
            bad_chars += 1
    print(f"[{label}] n={len(rows)}  length min/median/max = "
          f"{min(lens)}/{int(np.median(lens))}/{max(lens)}  (expected {expected_len})")
    print(f"[{label}] wrong length: {bad_len}   non-ACGTU chars: {bad_chars}   lowercase: {lower}")
    return bad_len == 0 and bad_chars == 0


def to_rna(seq):
    return seq.upper().replace("T", "U")


def fold_mfe(seq):
    structure, mfe = RNA.fold(to_rna(seq))
    return mfe


def dinuc_shuffle(seq):
    """Simple mononucleotide shuffle (keeps composition)."""
    s = list(seq)
    random.shuffle(s)
    return "".join(s)


def random_seq(length, gc=0.41):
    out = []
    for _ in range(length):
        if random.random() < gc:
            out.append(random.choice("GC"))
        else:
            out.append(random.choice("AU"))
    return "".join(out)


def summarize(name, E, length):
    E = np.asarray(E, dtype=float)
    print(f"  {name:28s} n={len(E):5d}  mean={E.mean():9.2f}  median={np.median(E):9.2f}  "
          f"min={E.min():9.2f}  max={E.max():9.2f}  per-nt mean={E.mean()/length:7.3f}  "
          f"#>0: {(E > 0).sum()}")


def cli_fold(seq):
    """Fold with the RNAfold command line tool and parse the MFE from the 2nd line."""
    p = subprocess.run(["RNAfold", "--noPS"], input=to_rna(seq) + "\n",
                       capture_output=True, text=True, check=True)
    lines = p.stdout.strip().split("\n")
    # lines[0] = sequence, lines[1] = "structure (-12.30)"
    tail = lines[1].rsplit("(", 1)[1].rstrip(")").strip()
    return float(tail)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("files", nargs="+", help="window files (chrom pos score seq)")
    ap.add_argument("--size", type=int, required=True, help="expected window length")
    ap.add_argument("--n", type=int, default=200, help="max sequences per file to fold")
    ap.add_argument("--cli", action="store_true", help="cross-check a few with RNAfold binary")
    ap.add_argument("--energies", default=None,
                    help="your own energies for the first file, one per line")
    args = ap.parse_args()

    random.seed(0)
    L = args.size
    all_ok = True

    for path in args.files:
        print("=" * 80)
        print(f"FILE: {path}")
        rows = read_windows(path, args.n)
        if not rows:
            print("  no rows read")
            all_ok = False
            continue

        ok = check_seqs(rows, L, path)
        all_ok &= ok

        seqs = [r[3] for r in rows]
        gc = np.mean([(s.upper().count("G") + s.upper().count("C")) / max(len(s), 1) for s in seqs])
        print(f"  mean GC content: {gc:.3f}")

        E_real = [fold_mfe(s) for s in seqs]
        E_shuf = [fold_mfe(dinuc_shuffle(s)) for s in seqs]
        E_rand = [fold_mfe(random_seq(L)) for _ in range(len(seqs))]

        print("\n  MFE summary (kcal/mol):")
        summarize("real windows", E_real, L)
        summarize("shuffled (same composition)", E_shuf, L)
        summarize(f"random GC=0.41 length {L}", E_rand, L)

        if any(e > 0 for e in E_real):
            print("  [FAIL] positive MFE found. MFE must be <= 0, so this is a parsing/computation bug.")
            all_ok = False
        if np.mean(E_real) / L < -1.0:
            print("  [WARN] mean MFE per nt < -1.0. That is unrealistically strong for length "
                  f"{L}; check whether sequences are longer than expected.")
            all_ok = False

        # Expected scaling check: MFE should be roughly proportional to length
        expected_mean = np.mean(E_rand)
        print(f"\n  Random baseline for length {L}: {expected_mean:.1f} kcal/mol")

        # CLI cross-check
        if args.cli:
            print("\n  CLI vs python API (first 5):")
            for s in seqs[:5]:
                try:
                    print(f"    python={fold_mfe(s):8.2f}   RNAfold={cli_fold(s):8.2f}")
                except Exception as e:
                    print(f"    CLI failed: {e}")
                    break

        # Compare to user's energies
        if args.energies and path == args.files[0]:
            mine = [float(x) for x in open(args.energies).read().split() if x.strip()]
            m = min(len(mine), len(E_real))
            diff = np.array(mine[:m]) - np.array(E_real[:m])
            print(f"\n  Your energies vs recomputed (first {m}):")
            print(f"    max |diff| = {np.abs(diff).max():.2f}   mean diff = {diff.mean():.2f}")
            if np.abs(diff).max() > 0.5:
                print("    [FAIL] your pipeline disagrees with a direct fold. Inspect rows with biggest diff:")
                for i in np.argsort(-np.abs(diff))[:5]:
                    print(f"      row {i}: yours={mine[i]:.2f} recomputed={E_real[i]:.2f} len={len(seqs[i])}")
                all_ok = False

        # Per-row preview
        print("\n  First 5 rows:")
        for (chrom, pos, score, seq), e in list(zip(rows, E_real))[:5]:
            print(f"    {chrom}:{pos} score={score} len={len(seq)} MFE={e:.2f}")

    print("=" * 80)
    print("ALL CHECKS PASSED" if all_ok else "SOME CHECKS FAILED, see messages above")


if __name__ == "__main__":
    main()
