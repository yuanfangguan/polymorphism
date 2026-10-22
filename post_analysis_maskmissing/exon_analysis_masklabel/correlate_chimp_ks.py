import glob
import numpy as np

record_ks={}
record_ka={}
FILE=open('human_chim.txt','r')
for line in FILE:
            line=line.strip()
            table=line.split('\t')
            try:
                record_ks[table[5]]=float(table[9]) ## ks value
                record_ka[table[5]]=float(table[10]) ## ks value
            except:
                pass
FILE.close()

FILE=open('gene_consolidated_mask.txt','r')
ks_vect=[]
ka_vect=[]
pred_vect=[]
for line in FILE:
    line=line.strip()
    table=line.split('\t')
    if (table[0] in record_ks):
        pred_vect.append(float(table[2]))
        ks_vect.append(record_ks[table[0]])
        ka_vect.append(record_ka[table[0]])
        print(record_ks[table[0]],table[2])
        
FILE.close()
ks_vect=np.asarray(ks_vect)
ka_vect=np.asarray(ka_vect)
pred_vect=np.asarray(pred_vect)
cor=np.corrcoef((ks_vect, pred_vect))
print(cor)
cor=np.corrcoef((ka_vect, pred_vect))
print(cor)


from scipy import stats
#res = stats.spearmanr(ks_vect, pred_vect, alternative='greater')
res = stats.spearmanr(ks_vect, pred_vect)#, alternative='greater')
print(res)
res = stats.spearmanr(ka_vect, pred_vect)#, alternative='greater')
print(res)
res = stats.spearmanr(ka_vect, ks_vect)#, alternative='greater')
print(res)
