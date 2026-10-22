import numpy as np
from sklearn import metrics
import pickle
import random
import tensorflow as tf
import sys
#from keras.utils.generic_utils import get_custom_objects
from keras.models import load_model

random.seed(3.14)

## same as cross entropy
def cross(y_true, y_pred):
    y_true = K.flatten(y_true)
    y_pred = tf.clip_by_value(K.flatten(y_pred), 1e-7, (1.0 - 1e-7))
    out = -(y_true * K.log(y_pred) + (1.0 - y_true) * K.log(1.0 - y_pred))
    return K.mean(out)


size=2048
channel=4

# first step 
model0=load_model('model0',custom_objects={'self_crossentropy': cross})
model1=load_model('model1',custom_objects={'self_crossentropy': cross})
model2=load_model('model2',custom_objects={'self_crossentropy': cross})
model3=load_model('model3',custom_objects={'self_crossentropy': cross})
model4=load_model('model4',custom_objects={'self_crossentropy': cross})

def load_chromosomes(filename):
    with open(filename, 'r') as file:
        chromosomes = [line.strip() for line in file.readlines()]
    return chromosomes

chr_list= load_chromosomes('test.chr')
chr_list_length=len(chr_list)


input_data={}
input_label={}
total_len_train={}
for the_chr in chr_list:
    input_data[the_chr]=pickle.load(open('../../data/'+the_chr+'.pkl', "rb"))
    input_label[the_chr]=pickle.load(open('../../data/'+the_chr+'_label_nonindel.pkl', "rb"))
    total_len_train[the_chr]=input_data[the_chr].shape[0]-2048


i=0
all_label=[]
all_sample=[]
while (i<50000):
#while (i<10):
    chr_id=int(random.random()*chr_list_length)
    the_chr=chr_list[chr_id]
    start=int(random.random()*total_len_train[the_chr])
    sample=input_data[the_chr][start:start+2048,:] 
    label=input_label[the_chr][start:start+2048] # minimal 0, maximal 2048
    if (np.min(np.max(sample,axis=1))==0):
        pass
    elif (np.sum(label)==0):
            pass
    else:
        all_sample.append(sample)
        all_label.append(label)
        i=i+1

pred0=np.asarray(model0.predict(np.asarray(all_sample)))
pred1=np.asarray(model1.predict(np.asarray(all_sample)))
pred2=np.asarray(model2.predict(np.asarray(all_sample)))
pred3=np.asarray(model3.predict(np.asarray(all_sample)))
pred4=np.asarray(model4.predict(np.asarray(all_sample)))
print(pred0.shape)

pred0=pred0.reshape(pred0.shape[0],pred0.shape[1])
pred1=pred1.reshape(pred1.shape[0],pred1.shape[1])
pred2=pred2.reshape(pred2.shape[0],pred2.shape[1])
pred3=pred3.reshape(pred3.shape[0],pred3.shape[1])
pred4=pred4.reshape(pred4.shape[0],pred4.shape[1])
pred=pred0+pred1+pred2+pred3+pred4
pred=pred/5 # n*2048

all_label=np.asarray(all_label)

print(all_label.shape)

all_auc=[]
all_cor=[]
i=0
while (i<2048):
    pred_tmp=pred[:,i]
    label_tmp=all_label[:,i]
    valid_indices = ~np.isnan(pred_tmp) & ~np.isnan(label_tmp)
    pred_tmp = pred_tmp[valid_indices]
    label_tmp = label_tmp[valid_indices]

    fpr, tpr, thresholds = metrics.roc_curve(label_tmp.flatten(),pred_tmp.flatten(), pos_label=1)
    all_auc.append(metrics.auc(fpr, tpr))
    cor=np.corrcoef(label_tmp,pred_tmp)
    all_cor.append(cor[0,1])

    i=i+1

np.savetxt('cor.perposition',np.asarray(all_cor))
np.savetxt('auc.perposition',np.asarray(all_auc))
