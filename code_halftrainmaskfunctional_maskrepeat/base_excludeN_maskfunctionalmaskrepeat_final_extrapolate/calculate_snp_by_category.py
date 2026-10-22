import numpy as np
import pickle
import glob

## record snp prediction
snp={}
all_bin=glob.glob('NC*.npy')
for the_label in all_bin:
    data= np.load(the_label)
    name=the_label.split('.pred.npy')[0]#.split('/')[-1]
    snp[name]=data

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



TYPE=open('all_types.txt','r')
type_map={}
type_i=0
for line in TYPE:
    line=line.strip()
    type_map[line]=type_i
    type_i=type_i+1
TYPE.close()

sum_length={}
sum_ones={}

GENOME=open('/local/disk4/gyuanfan/references/GRCh38_latest_genomic.gff','r')
for line in GENOME:
    if ('#' in line):
        pass
    else:
        line=line.strip()
        table=line.split('\t')
        if (table[1] == 'RefSeq'):
            pass
        else:
                
#            print(table[0],table[3],table[4])
            if (table[3] == table[4]):
                pass

            elif (table[0] in map_genome):
                if (np.min(map_genome[table[0]][int(table[3]):int(table[4])])<0):
                    print('excluded')
                    pass
                else:

                    if ((table[0] in snp) and (table[2] in type_map)):
                        if (table[2] in sum_length):
                            sum_length[table[2]]=sum_length[table[2]]+(int(table[4])-int(table[3]))
                            sum_ones[table[2]]=sum_ones[table[2]]+np.sum(snp[table[0]][int(table[3]):int(table[4])])
                        else:
                            sum_length[table[2]]=int(table[4])-int(table[3])
                            sum_ones[table[2]]=np.sum(snp[table[0]][int(table[3]):int(table[4])])
GENOME.close()

all_types=type_map.keys()

NEW=open('snp_by_type_bin.txt','w')
for the_type in all_types:
    NEW.write(the_type)
    NEW.write('\t')
    try:
        avg=sum_ones[the_type]/sum_length[the_type]
    except:
        avg='NA'
    NEW.write(str(avg))
    NEW.write('\n')




