import numpy as np
import pickle
import sys



map_ACTG={}
map_ACTG['A']=[1,0]
map_ACTG['C']=[0,1]
map_ACTG['T']=[0,0]
map_ACTG['G']=[1,1]


GENOME=open('GRCh38_latest_genomic.fna','r')
#map_genome=np.zeros((250000000,4))
map_genome=np.zeros((25000000,2)) # one tenth of genome
i=0
for line in GENOME:
    line=line.strip()
    if ('>' in line):
        table=line.split(' ')
        table=table[0].split('>')
        the_chr=table[1]
    else:
        line=line.upper()
        for aaa in line:
            if (the_chr ==sys.argv[1]):
                try:
                    map_genome[i,:]=map_ACTG[aaa]
                except:
                    pass
                i=i+1
GENOME.close()



print(the_chr)
pickle.dump(map_genome, open(sys.argv[1]+'.bit.pkl', 'wb'))
