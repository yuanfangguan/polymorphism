import numpy as np
import pickle
import glob
import os

## Define chromosome mapping
map_chr = {
    'NC_000001.11': 'chr1', 'NC_000002.12': 'chr2', 'NC_000003.12': 'chr3', 
    'NC_000004.12': 'chr4', 'NC_000005.10': 'chr5', 'NC_000006.12': 'chr6', 
    'NC_000007.14': 'chr7', 'NC_000008.11': 'chr8', 'NC_000009.12': 'chr9', 
    'NC_000010.11': 'chr10', 'NC_000011.10': 'chr11', 'NC_000012.12': 'chr12', 
    'NC_000013.11': 'chr13', 'NC_000014.9': 'chr14', 'NC_000015.10': 'chr15', 
    'NC_000016.10': 'chr16', 'NC_000017.11': 'chr17', 'NC_000018.10': 'chr18', 
    'NC_000019.10': 'chr19', 'NC_000020.11': 'chr20', 'NC_000021.9': 'chr21', 
    'NC_000022.11': 'chr22', 'NC_000023.11': 'chrX', 'NC_000024.10': 'chrY'
}
map_chr_back = {v: k for k, v in map_chr.items()}

## Read genome sequence
GENOME = open('/local/disk4/gyuanfan/references/GRCh38_latest_genomic.fna', 'r')
map_ACTG = {'A': [1, 0, 0, 0], 'C': [0, 1, 0, 0], 'T': [0, 0, 1, 0], 'G': [0, 0, 0, 1]}
map_genome = {}

for line in GENOME:
    line = line.strip()
    if '>' in line:
        table = line.split(' ')[0].split('>')
        the_chr = table[1]
        if 'NC_' in the_chr:
            map_genome[the_chr] = []
    else:
        line = line.upper()
        if 'NC_' in the_chr:
            map_genome[the_chr].extend(1 if base in map_ACTG else -1 for base in line)

GENOME.close()

## Convert lists to numpy arrays
for the_chr in map_genome.keys():
    map_genome[the_chr] = np.asarray(map_genome[the_chr])

## Load SNP data
snp = {}
all_bin=glob.glob('../../calibrated_predictions_masklabel/*.pred.npy')

for the_label in all_bin:
    if 'non' not in the_label:
        with open(the_label, 'rb') as pickle_file:
            data = np.load(pickle_file)
        name = the_label.split('.pred.npy')[0].split('/')[-1]
        print(name)
        snp[name] = data

## Read types
TYPE = open('all_types.txt', 'r')
type_map = {line.strip(): idx for idx, line in enumerate(TYPE)}
TYPE.close()

sum_length = {}
sum_ones = {}

## Create output directory
output_dir = "predicted_output_segments"
os.makedirs(output_dir, exist_ok=True)

## Initialize a dictionary to store SNP values by type
feature_type_data = {feature: [] for feature in type_map.keys()}

## Process genomic annotations
GENOME = open('/local/disk4/gyuanfan/references/GRCh38_latest_genomic.gff', 'r')

for line in GENOME:
    if not line.startswith('#'):
        table = line.strip().split('\t')
        if table[1] != 'RefSeq' and table[0] in snp and table[2] in type_map:
            start, end = int(table[3]), int(table[4])
            if start != end and np.min(map_genome[table[0]][start:end]) >= 0:
                segment_value = np.sum(snp[table[0]][start:end])

                # Store segment values by annotation type
                feature_type_data[table[2]].append(segment_value/(end-start))
                # Aggregate statistics
                sum_length[table[2]] = sum_length.get(table[2], 0) + (end - start)
                sum_ones[table[2]] = sum_ones.get(table[2], 0) + segment_value

GENOME.close()

## Save SNP segment values by annotation type
for feature_type, segments in feature_type_data.items():
    print(feature_type, len(segments))
    if segments:  # Only save if there are any segments for this feature type
        type_filename = os.path.join(output_dir, f"{feature_type}.txt")
        with open(type_filename, 'w') as type_file:
            for segment in segments:
                type_file.write(f"{segment}\n")

## Write overall SNP type averages
with open('snp_by_type_bin.txt', 'w') as NEW:
    for the_type in type_map.keys():
        avg = sum_ones.get(the_type, 0) / sum_length[the_type] if sum_length.get(the_type, 0) > 2000 else 'NA'
        NEW.write(f"{the_type}\t{avg}\n")

