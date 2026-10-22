import pyBigWig

bw = pyBigWig.open("RNAseq_default.bw")
print(bw.chroms())  # Print available chromosome names and their sizes
bw.close()

