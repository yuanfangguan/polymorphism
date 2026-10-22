import numpy as np
import pickle
import sys


## set label length
map_chr={}
map_chr['NC_000001.11']='chr1'
map_chr['NC_000002.12']='chr2'
map_chr['NC_000003.12']='chr3'
map_chr['NC_000004.12']='chr4'
map_chr['NC_000005.10']='chr5'
map_chr['NC_000006.12']='chr6'
map_chr['NC_000007.14']='chr7'
map_chr['NC_000008.11']='chr8'
map_chr['NC_000009.12']='chr9'
map_chr['NC_000009.12']='chr9'
map_chr['NC_000010.11']='chr10'
map_chr['NC_000011.10']='chr11'
map_chr['NC_000012.12']='chr12'
map_chr['NC_000013.11']='chr13'
map_chr['NC_000014.9']='chr14'
map_chr['NC_000015.10']='chr15'
map_chr['NC_000016.10']='chr16'
map_chr['NC_000017.11']='chr17'
map_chr['NC_000018.10']='chr18'
map_chr['NC_000019.10']='chr19'
map_chr['NC_000020.11']='chr20'
map_chr['NC_000021.9']='chr21'
map_chr['NC_000022.11']='chr22'
map_chr['NC_000023.11']='chrX'
map_chr['NC_000024.10']='chrY'


all_chr=map_chr.keys()
map_chr_back={}
for the_chr in all_chr:
    map_chr_back[map_chr[the_chr]]=the_chr



GENOME=open('GRCh38_latest_genomic.fna','r')
map_ACTG={}
map_ACTG['A']=[1,0,0,0]
map_ACTG['C']=[0,1,0,0]
map_ACTG['T']=[0,0,1,0]
map_ACTG['G']=[0,0,0,1]
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

all_chr=map_genome.keys()
label={}
for the_chr in all_chr:
    map_genome[the_chr]=np.asarray(map_genome[the_chr]).shape[0]
    #print(map_genome[the_chr])
    label[the_chr]=np.zeros(map_genome[the_chr])

LABEL=open('1000genome/1000_genome_snv.dat','r')
for line in LABEL:
    line=line.strip()
    table=line.split('\t')
    label[map_chr_back[table[0]]][int(table[1])-1:int(table[2])-1]=1 # whether there is a mutation
LABEL.close()

for the_chr in all_chr:
    pickle.dump(label[the_chr], open(the_chr+'_label_nonindel.pkl', 'wb'))
