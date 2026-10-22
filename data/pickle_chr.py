import numpy as np
import pickle
import sys

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
for the_chr in all_chr:
        map_genome[the_chr]=np.asarray(map_genome[the_chr])
        pickle.dump(map_genome[the_chr], open(the_chr+'.pkl', 'wb'), protocol=4)
