import random
import sys
random.seed(sys.argv[1])
FILE=open('/local/disk4/gyuanfan/references/GRCh38_latest_genomic.gff','r')
the_type={}
for line in FILE:
    table=line.split('\t')
    try:
        if ('biological_region' ==table[2]):
            pass
        elif ('region' == table[2]):
            pass
        elif ('chromosome' == table[2]):
            pass
        elif ('conserved_region'==table[2]):
            pass
        elif ('sequence_comparison' == table[2]):
            pass
        elif ('sequence_feature'==table[2]):
            pass
        elif ('sequence_alteration' == table[2]):
            pass
        elif ('match' == table[2]):
            pass
        elif ('scaffold' == table[2]):
            pass
        else:
            the_type[table[2]]=1
    except:
        pass
FILE.close()

all_type=the_type.keys()
TRAIN=open('train_type.txt','w')
TEST=open('test_type.txt','w')

for the_type in all_type:
    rrr=random.random()
    if (rrr<0.5):
        TRAIN.write(the_type)
        TRAIN.write('\n')
    else:
        TEST.write(the_type)
        TEST.write('\n')
TRAIN.close()
TEST.close()

GENOME=open('/local/disk4/gyuanfan/references/GRCh38_latest_genomic.fna','r')
map_genome={}
for line in GENOME:
    line=line.strip()
    if ('>' in line):
        if ('unlocalized' in line):
            pass
        else:
            if ('unplaced' in line):
                pass
            else:

                if ('alternate' in line):
                    pass
                else:
                    if ('mitochondrion' in line):
                        pass
                    else:
                        if ('patch' in line):
                            pass
                        else:
                            if ('chromosome X' in line) or ('chromosome Y' in line):
                                pass
                            else:
                                print(line)
                                table=line.split(' ')
                                table=table[0].split('>')
                                the_chr=table[1]
                                map_genome[the_chr]=[]
GENOME.close()

all_chr=map_genome.keys()
TRAIN=open('train_chr.txt','w')
TEST=open('test_chr.txt','w')
i=0
for the_chr in all_chr:
#    if (i%10)==0:
        TRAIN.write(the_chr)
        TRAIN.write('\n')
        TEST.write(the_chr)
        TEST.write('\n')
TRAIN.close()
TEST.close()

