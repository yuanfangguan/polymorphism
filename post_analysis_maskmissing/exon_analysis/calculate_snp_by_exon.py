import numpy as np
### record snp binary number
import pickle
import glob

## record N areas

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



GENOME=open('/local/disk4/gyuanfan/references/GRCh38_latest_genomic.fna','r')
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
                  if (aaa in map_ACTG):
                      map_genome[the_chr].append(1)
                  
                  else:
                      map_genome[the_chr].append(-1)
GENOME.close()


all_chr=map_genome.keys()
for the_chr in all_chr:
    map_genome[the_chr]=np.asarray(map_genome[the_chr])


snp={}
all_bin=glob.glob('/local/disk4/gyuanfan/rate/data/*_label.pkl')
for the_label in all_bin:
    if ('non' in the_label):
        pass
    else:
        with open(the_label, 'rb') as pickle_file:
            data= pickle.load(pickle_file)
        name=the_label.split('_label.pkl')[0].split('/')[-1]
        snp[name]=data



GENOME=open('/local/disk4/gyuanfan/references/GRCh38_latest_genomic.gff','r')
NEW=open('exon_observed.txt','w')
for line in GENOME:
    if ('#' in line):
        pass
    else:
        line=line.strip()
        table=line.split('\t')
        if ('RefSeq' in table[1]):
            if ((table[0] in snp) and (table[2] =='exon')):
                if ('pseudogene' in line):
                    pass
                else:
                    if (table[3] == table[4]):
                        pass
                    elif (np.min(map_genome[table[0]][int(table[3]):int(table[4])]<0)):
                        pass
                    else:
                        NEW.write(line)
                        NEW.write('\t')
                        ratio=np.sum(snp[table[0]][int(table[3]):int(table[4])])/(int(table[4])-int(table[3]))
                        NEW.write(str(ratio))
                        NEW.write('\n')
GENOME.close()

