import glob

all_files=glob.glob("*vcf")
INDEL=open('1000_genome_indel.dat','w')
SNV=open('1000_genome_snv.dat','w')
BOTH=open('1000_genome_all.dat','w')


for the_file in all_files:
    FILE=open(the_file,'r')
    for line in FILE:
        if '#' in line:
            pass
        else:
            table=line.split('\t')
            if (len(table[4]) == len(table[3])):
                    SNV.write(table[0]+'\t'+table[1]+'\t'+str(int(table[1])+1)+'\n')
                    BOTH.write(table[0]+'\t'+table[1]+'\t'+str(int(table[1])+1)+'\n')
            else:
                    INDEL.write(table[0]+'\t'+table[1]+'\t'+str(int(table[1])+len(table[3])+1)+'\n')
                    BOTH.write(table[0]+'\t'+table[1]+'\t'+str(int(table[1])+len(table[3])+1)+'\n')
    FILE.close()

