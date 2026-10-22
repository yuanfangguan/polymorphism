import os
import numpy as np
from Bio import SeqIO
import heapq  # 用于高效维护 top 100
import pyBigWig

# --- 路径配置 ---
genome_fasta = "../../../../references/GRCh38_latest_genomic.fna"
gff_file = "../../../../references/GRCh38_latest_genomic.gff"
phylop_bw_file = "hg38.phyloP100way.bw"  # 确保你已经下载了这个文件
window_size = 64
step_size = 32
top_k = 1000

# NCBI 到 UCSC 染色体名字转换
def ncbi_to_ucsc(ncbi_id):
    try:
        num = int(ncbi_id.split('.')[0].split('_')[1])
        if num == 23: return "chrX"
        if num == 24: return "chrY"
        if 1 <= num <= 22: return f"chr{num}"
    except:
        pass
    return None

def get_highest_phylop_windows(genome_fasta, map_genome, bw_file):
    # Python 默认是小顶堆 (min-heap)。
    # 我们把 PhyloP 得分存进去，堆顶永远是当前 100 个里最小的。
    # 如果新窗口的分数比堆顶大，就把堆顶弹出来，把新窗口加进去，这样最后留下的就是最大的 100 个。
    highest_windows = []
    bw = pyBigWig.open(bw_file)

    for record in SeqIO.parse(genome_fasta, "fasta"):
        chrom_name = record.id
        print(f"Processing {chrom_name}...")

        ucsc_chrom = ncbi_to_ucsc(chrom_name)
        if not ucsc_chrom or ucsc_chrom not in bw.chroms():
            print(f"Chromosome {chrom_name} ({ucsc_chrom}) not found in PhyloP bigWig. Skipping.")
            continue

        if chrom_name not in map_genome:
            continue

        sequence = str(record.seq).upper()

        for i in range(0, len(sequence) - window_size + 1, step_size):
            # 1. 检查是否在 mask 区域 (这部分和你原代码逻辑完全一致)
            mask_window = map_genome[chrom_name][i:i + window_size, :]
            if np.any(np.all(mask_window == [0, 0, 0, 0], axis=1)):
                continue

            window_seq = sequence[i:i + window_size]
            # 忽略含有 N 的模糊序列
            if 'N' in window_seq:
                continue

            # 2. 从 BigWig 提取该窗口的平均 PhyloP 得分
            try:
                # 提取均值
                stats = bw.stats(ucsc_chrom, i, i + window_size, type="mean")
                if stats[0] is None:
                    continue
                avg_phylop = float(stats[0])
            except RuntimeError:
                continue

            # 3. 存入小顶堆维护 Top 100 最高分
            if len(highest_windows) < top_k:
                heapq.heappush(highest_windows, (avg_phylop, chrom_name, i, window_seq))
            else:
                # 如果当前得分大于堆中最小的那个，就替换掉它
                if avg_phylop > highest_windows[0][0]:
                    heapq.heappushpop(highest_windows, (avg_phylop, chrom_name, i, window_seq))

    bw.close()
    return highest_windows

if __name__ == "__main__":
    # ---------------- 你的原始 Mask 构建逻辑开始 ----------------
    print("Building mask array...")
    GENOME=open(genome_fasta,'r')
    map_ACTG = {'A':[1,0,0,0], 'C':[0,1,0,0], 'T':[0,0,1,0], 'G':[0,0,0,1]}
    map_genome={}
    
    for line in GENOME:
        line=line.strip()
        if ('>' in line):
            table=line.split(' ')
            table=table[0].split('>')
            the_chr=table[1]
            if ('NC_' in the_chr):
                map_genome[the_chr]=[]
        else:
            line=line.upper()
            for aaa in line:
                if ('NC_' in the_chr):
                    try:
                        map_genome[the_chr].append(map_ACTG[aaa])
                    except:
                        map_genome[the_chr].append([0,0,0,0])
    GENOME.close()

    all_chr = list(map_genome.keys())
    for the_chr in all_chr:
        map_genome[the_chr]=np.asarray(map_genome[the_chr])

    print("Applying GFF masks...")
    GFF=open(gff_file,'r')
    for line in GFF:
        if ('#' in line):
            pass
        else:
            line=line.strip()
            table=line.split('\t')
            if (table[1] == 'RefSeq'):
                pass
            else:
                if (table[0] in map_genome):
                    try:
                        # 你的保守做法：注释区域前后 512 碱基全部排查
                        map_genome[table[0]][(int(table[3])-1)-512:int(table[4])+512,:]=0
                    except:
                        pass
    GFF.close()
    print("Masking complete.")
    # ---------------- 你的原始 Mask 构建逻辑结束 ----------------

    # 开始寻找最高 PhyloP 的 100 个窗口
    print("Scanning genome for highest PhyloP scores...")
    highest_windows = get_highest_phylop_windows(genome_fasta, map_genome, phylop_bw_file)

    # 按分数从高到低排序
    highest_windows = sorted(highest_windows, key=lambda x: x[0], reverse=True)

    output_file = "highest_100_64_phylop_windows.txt"
    with open(output_file, "w") as f:
        for avg_phylop, chrom, pos, seq in highest_windows:
            f.write(f"{chrom}\t{pos}\t{avg_phylop:.4f}\t{seq}\n")

    print(f"Top 100 highest PhyloP windows saved to {output_file}")
