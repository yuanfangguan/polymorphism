

REF=open('observed_predict_comparison_mask.txt','r')
length={}
count_obs={}
count_pred={}
for line in REF:
    line=line.strip()
    table=line.split('\t')
    ttt=table[3].split('gene=')
    ttt=ttt[1].split(';')
    if ttt[0] in length:
        length[ttt[0]]=length[ttt[0]]+(int(table[2])-int(table[1]))
        count_obs[ttt[0]]=count_obs[ttt[0]]+(int(table[2])-int(table[1]))*float(table[4])
        count_pred[ttt[0]]=count_pred[ttt[0]]+(int(table[2])-int(table[1]))*float(table[5])
    else:
        length[ttt[0]]=(int(table[2])-int(table[1]))
        count_obs[ttt[0]]=(int(table[2])-int(table[1]))*float(table[4])
        count_pred[ttt[0]]=(int(table[2])-int(table[1]))*float(table[5])

REF.close()

all_genes=length.keys()
NEW=open('gene_consolidated_mask.txt','w')
NEW.write('gene name\tobserved\tpredicted\n')
for the_gene in all_genes:
    NEW.write(the_gene)
    NEW.write('\t')

    ooo1=count_obs[the_gene]/length[the_gene]
    NEW.write(str(ooo1))
    NEW.write('\t')
    ooo2=count_pred[the_gene]/length[the_gene]
    NEW.write(str(ooo2))
    NEW.write('\t')
    try:
        ratio=ooo1/ooo2
        NEW.write(str(ratio))
    except:
        print(the_gene)
        pass

    NEW.write('\n')
NEW.close()
