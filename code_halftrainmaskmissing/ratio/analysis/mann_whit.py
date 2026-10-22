from scipy.stats import mannwhitneyu
import sys

FILE1=open(sys.argv[1],'r')
FILE2=open(sys.argv[2],'r')
group1=[]
group2=[]
for line in FILE1:
    val=float(line.strip())
    group1.append(val)
FILE1.close()
for line in FILE2:
    val=float(line.strip())
    group2.append(val)
FILE2.close()
# Example data
#group1 = [12, 15, 14, 10, 9, 8, 16]  # Sample 1
#group2 = [22, 25, 20, 18, 19, 23, 21]  # Sample 2


# Perform Mann-Whitney U test
stat, p_value = mannwhitneyu(group1, group2, alternative='two-sided')

print(f"Mann-Whitney U statistic: {stat}")
print(f"P-value: {p_value}")

