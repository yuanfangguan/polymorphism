

REF=open('exon_predicted_mask.txt','r')
predict={}
for line in REF:
    line=line.strip()
    table=line.split('\t')
    concat=table[0]+'_'+table[3]+'_'+table[4]
    predict[concat]=table[9]

    
REF.close()

FILE=open('exon_observed.txt','r')
NEW=open('observed_predict_comparison_mask.txt','w')
for line in FILE:
    line=line.strip()
    table=line.split('\t')
    concat=table[0]+'_'+table[3]+'_'+table[4]
    if (concat in predict):
        NEW.write(table[0])
        NEW.write('\t')
        NEW.write(table[3])
        NEW.write('\t')
        NEW.write(table[4])
        NEW.write('\t')
        NEW.write(table[8])
        NEW.write('\t')
        NEW.write(table[9])
        NEW.write('\t')

        NEW.write(predict[concat])
        NEW.write('\n')
NEW.close()
FILE.close()


