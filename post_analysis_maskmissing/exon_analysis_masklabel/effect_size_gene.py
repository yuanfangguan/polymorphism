import numpy as np

observed=[]
predicted=[]
FILE=open('gene_consolidated.txt','r')
line=FILE.readline()
for line in FILE:
    line=line.strip()
    table=line.split('\t')
    if (float(table[1])==0) or (float(table[2])==0):
        pass
    else:
        observed.append(float(table[1]))
        predicted.append(float(table[2]))
FILE.close()
# Sorting observed and predicted values based on observed values
observed_sorted, predicted_sorted = zip(*sorted(zip(observed, predicted)))

# Finding the index corresponding to the 95th percentile
index_95th_percentile = int(0.95 * len(observed_sorted))

# Extracting the top 5% observed values and their corresponding predicted values
top_5_percent_observed = observed_sorted[index_95th_percentile:]
top_5_percent_predicted = predicted_sorted[index_95th_percentile:]

# Calculating median of top 5% observed values
if len(top_5_percent_observed) % 2 == 0:
    median_observed = (top_5_percent_observed[len(top_5_percent_observed) // 2] +
                       top_5_percent_observed[len(top_5_percent_observed) // 2 - 1]) / 2
else:
    median_observed = top_5_percent_observed[len(top_5_percent_observed) // 2]

# Calculating median of corresponding predicted values
if len(top_5_percent_predicted) % 2 == 0:
    median_predicted = (top_5_percent_predicted[len(top_5_percent_predicted) // 2] +
                         top_5_percent_predicted[len(top_5_percent_predicted) // 2 - 1]) / 2
else:
    median_predicted = top_5_percent_predicted[len(top_5_percent_predicted) // 2]

print("Median value of top 5% observed:", median_observed)
print("Median value of corresponding predicted values:", median_predicted)


# Sorting observed and predicted values based on observed values
observed_sorted, predicted_sorted = zip(*sorted(zip(observed, predicted)))
observed = []
predicted = []

# Open the file and read data into lists
with open('gene_consolidated.txt', 'r') as file:
    line=file.readline()
    for line in file:
        line = line.strip()
        table = line.split('\t')
        if (float(table[1])==0) or (float(table[2])==0):
            pass
        else:
            observed.append(float(table[1]))
            predicted.append(float(table[2]))

# Sorting observed and predicted values based on observed values
observed_sorted, predicted_sorted = zip(*sorted(zip(observed, predicted)))

# Finding the index corresponding to the 5th percentile (bottom 5%)
index_5th_percentile = int(0.05 * len(observed_sorted))

# Extracting the bottom 5% observed values and their corresponding predicted values
bottom_5_percent_observed = observed_sorted[:index_5th_percentile]
bottom_5_percent_predicted = predicted_sorted[:index_5th_percentile]

# Calculating median of bottom 5% observed values
if len(bottom_5_percent_observed) % 2 == 0:
    median_bottom_observed = (bottom_5_percent_observed[len(bottom_5_percent_observed) // 2] +
                              bottom_5_percent_observed[len(bottom_5_percent_observed) // 2 - 1]) / 2
else:
    median_bottom_observed = bottom_5_percent_observed[len(bottom_5_percent_observed) // 2]

# Calculating median of corresponding predicted values
if len(bottom_5_percent_predicted) % 2 == 0:
    median_bottom_predicted = (bottom_5_percent_predicted[len(bottom_5_percent_predicted) // 2] +
                               bottom_5_percent_predicted[len(bottom_5_percent_predicted) // 2 - 1]) / 2
else:
    median_bottom_predicted = bottom_5_percent_predicted[len(bottom_5_percent_predicted) // 2]

print("Median value of bottom 5% observed:", median_bottom_observed)
print("Median value of corresponding predicted values:", median_bottom_predicted)


