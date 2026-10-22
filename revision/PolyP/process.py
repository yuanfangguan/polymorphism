#!/usr/bin/env python
import os, argparse, numpy as np, pyBigWig
from numpy.lib.format import open_memmap

def main():
    ap = argparse.ArgumentParser(description="Convert bigWig to per-chromosome .npy arrays")
    ap.add_argument("bw", help="input bigWig, e.g. hg38.phyloP4way.bw")
    ap.add_argument("outdir", help="output folder for per-chromosome .npy files")
    ap.add_argument("--chunk", type=int, default=1_000_000, help="chunk size (bases) to stream per read")
    ap.add_argument("--dtype", choices=["float32","float16"], default="float32", help="dtype of saved arrays")
    ap.add_argument("--chroms", help="comma-separated list; default = standard chrs present in file")
    ap.add_argument("--fill", type=float, default=None, help="replace NaN with this value; default keeps NaN")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    bw = pyBigWig.open(args.bw)
    chroms_dict = bw.chroms()  # {'chr1': length, ...}

    if args.chroms:
        wanted = [c.strip() for c in args.chroms.split(",")]
    else:
        std = [f"chr{i}" for i in range(1,23)] + ["chrX","chrY","chrM","chrMT"]
        wanted = [c for c in std if c in chroms_dict]

    dtype = np.float32 if args.dtype == "float32" else np.float16

    for chrom in wanted:
        L = chroms_dict[chrom]
        path = os.path.join(args.outdir, f"{chrom}.npy")
        mm = open_memmap(path, mode="w+", dtype=dtype, shape=(L,))  # 生成真正的 .npy 文件（可 mmap 读取）
        start = 0
        while start < L:
            end = min(start + args.chunk, L)
            vals = np.array(bw.values(chrom, start, end, numpy=True), dtype=np.float64)  # 包含 NaN
            if args.fill is not None:
                vals = np.nan_to_num(vals, nan=args.fill)
            mm[start:end] = vals.astype(dtype, copy=False)
            start = end
        mm.flush(); del mm
        print(f"Wrote {chrom}: {L:,} bases -> {path}")
    bw.close()

if __name__ == "__main__":
    main()

