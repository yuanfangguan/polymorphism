# record two sets of chromsomes:
# 1. before: original chromosome seq
# 2. after SNP edit chromosome seq

## record original chromosome seq
GENOME=open('/local/disk3/gyuanfan/references/GRCh38_latest_genomic.fna','r')
genome={}
for line in GENOME:
    line=line.strip()
    if ('>' in line):
          table=line.split(' ')
          table=table[0].split('>')
          the_chr=table[1]
          if ('NC_' in the_chr):
              genome[the_chr]=""
    else:
          line=line.upper()
          if ('NC_' in the_chr):
              genome[the_chr]=genome[the_chr]+line
GENOME.close()

mutated_genome=genome

REF=open('/local/disk3/gyuanfan/rate/data/dbSnp153.bed','r')
for line in REF:
    line=line.strip()
    table=line.split('\t')
    if (len(table[4])==len(table[6])):
        mutated_genome[table[0]][int(table[1]):int(table[2])]=table[6]
REF.close()


from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio.Align import AlignInfo

def calculate_ka_ks(seq1, seq2):
    # Define a function to count synonymous and nonsynonymous mutations
    def count_mutations(codon1, codon2):
        # Initialize counters
        synonymous = 0
        nonsynonymous = 0
        # Ensure codons are of the same length
        min_len = min(len(codon1), len(codon2))
        # Loop through each position in the codons
        for i in range(min_len):
            if codon1[i] != codon2[i]:
                # If the codons differ
                if Seq(codon1[i:i+3]).translate() == Seq(codon2[i:i+3]).translate():
                    # If the mutations are synonymous
                    synonymous += 1
                else:
                    # If the mutations are nonsynonymous
                    nonsynonymous += 1
        return synonymous, nonsynonymous

    # Initialize counts
    total_synonymous = 0
    total_nonsynonymous = 0

    # Iterate over codons in the aligned sequences
    for i in range(0, len(seq1), 3):
        codon1 = seq1[i:i+3]
        codon2 = seq2[i:i+3]
        synonymous, nonsynonymous = count_mutations(codon1, codon2)
        total_synonymous += synonymous
        total_nonsynonymous += nonsynonymous

    # Calculate Ka and Ks
    try:
        ka = total_nonsynonymous / (len(seq1) // 3)
        ks = total_synonymous / (len(seq1) // 3)
        ka_ks_ratio = ka / ks
        return ka_ks_ratio,ka,ks
    except ZeroDivisionError:
        return None, ka,ks

### go through all exons:
GENOME=open('/local/disk3/gyuanfan/references/GRCh38_latest_genomic.gff','r')
NEW=open('exon_kaks.txt','w')
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
                    else:
                        NEW.write(line)
                        NEW.write('\t')
                        seq1=genome[table[0]][int(table[3]):int(table[4])]
                        seq2=mutated_genome[table[0]][int(table[3]):int(table[4])]
                        ka_ks_ratio,ka,ks = calculate_ka_ks(seq1, seq2)
                        NEW.write(str(ka_ks_ratio))
                        NEW.write('\t')
                        NEW.write(str(ka))
                        NEW.write('\t')
                        NEW.write(str(ks))
                        NEW.write('\n')
GENOME.close()

