import sys
import pyBigWig
import pandas as pd

# Load the BigWig file
bw = pyBigWig.open("RNAseq_default.bw")

# Read the regions file (chromosome, start position)
regions = pd.read_csv("../ratio_masklabel/analysis/lowest_100_"+sys.argv[1]+"_windows.txt", sep="\t", header=None, names=["chrom", "start","value","seq"])

# Define the window size
window_size = int(sys.argv[1])

# List to store results
results = []

# Iterate over each row
for _, row in regions.iterrows():
    chrom = row["chrom"].replace("NC_000001.11", "chr1")
    chrom = chrom.replace("NC_000002.12", "chr2")
    chrom = chrom.replace("NC_000003.12", "chr3")
    chrom = chrom.replace("NC_000004.12", "chr4")
    chrom = chrom.replace("NC_000005.10", "chr5")
    chrom = chrom.replace("NC_000006.12", "chr6")
    chrom = chrom.replace("NC_000007.14", "chr7")
    chrom = chrom.replace("NC_000008.11", "chr8")
    chrom = chrom.replace("NC_000009.12", "chr9")
    chrom = chrom.replace("NC_000010.11", "chr10")
    chrom = chrom.replace("NC_000011.10", "chr11")
    chrom = chrom.replace("NC_000012.12", "chr12")
    chrom = chrom.replace("NC_000013.11", "chr13")
    chrom = chrom.replace("NC_000014.9", "chr14")
    chrom = chrom.replace("NC_000015.10", "chr15")
    chrom = chrom.replace("NC_000016.10", "chr16")
    chrom = chrom.replace("NC_000017.11", "chr17")
    chrom = chrom.replace("NC_000018.10", "chr18")
    chrom = chrom.replace("NC_000019.10", "chr19")
    chrom = chrom.replace("NC_000020.11", "chr20")
    chrom = chrom.replace("NC_000021.9", "chr21")
    chrom = chrom.replace("NC_000022.11", "chr22")
    chrom = chrom.replace("NC_000023.11", "chrX")
    chrom = chrom.replace("NC_000024.10", "chrY")


    start = int(row["start"])
    end = start + window_size

    # Extract values from the BigWig file
    try:
        values = bw.values(chrom, start, end, numpy=True)

    # Compute the average while ignoring NaN values
        avg_expression = values[~pd.isna(values)].mean() if values is not None else "NA"
        

    # Store result
        results.append([chrom, start, end, avg_expression])
    except:
        print(chrom,start,end)

# Close the BigWig file
bw.close()

# Convert results to a DataFrame
output_df = pd.DataFrame(results, columns=["chrom", "start", "end", "avg_expression"])

# Save to a file
output_df.to_csv(sys.argv[1]+"average_expression_results.txt", sep="\t", index=False)

print("Average expression values saved to average_expression_results.txt")

