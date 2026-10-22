import numpy as np
import matplotlib.pyplot as plt

# 读取第一个文件
file1 = "1024_mean_expression.txt"
data1 = []
with open(file1, "r") as f:
    for line in f:
        value = line.strip()
        data1.append(float(value) if value != "nan" else 0)  # 处理 'nan' 为 0

# 读取第二个文件
file2 = "whole_genome/1024average_expression_results.txt"
data2 = []
with open(file2, "r") as f:
    next(f)  # 跳过表头
    for line in f:
        columns = line.strip().split("\t")
        value = columns[-1] if len(columns) > 3 else "0"  # 如果缺失，填充为 0
        data2.append(float(value) if value != "nan" else 0)

# 进行 log(x+1) 变换
data1 = np.log1p(np.array(data1))  
data2 = np.log1p(np.array(data2))

# 画直方图
plt.figure(figsize=(8, 5))
plt.hist(data1, bins=50, alpha=0.5, color="blue", label="1024_mean_expression (log)", density=True)
plt.hist(data2, bins=50, alpha=0.5, color="orange", label="whole_genome_expression (log)", density=True)

plt.xlabel("Log(Avg Expression + 1)")
plt.ylabel("Density")
plt.legend()
plt.title("Comparison of Expression Distributions (Log-Transformed Histogram)")
plt.grid(True)

# 保存图片
plt.savefig("expression_histogram_log_comparison.png", dpi=300)
plt.show()

