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


size=1024
channel=4

# first step 
model0=load_model('model0',custom_objects={'self_crossentropy': cross})
model1=load_model('model1',custom_objects={'self_crossentropy': cross})
model2=load_model('model2',custom_objects={'self_crossentropy': cross})
model3=load_model('model3',custom_objects={'self_crossentropy': cross})
model4=load_model('model4',custom_objects={'self_crossentropy': cross})

chr_list=['NC_000001.11','NC_000002.12','NC_000004.12','NC_000006.12','NC_000008.11','NC_000010.11','NC_000012.12','NC_000014.9','NC_000016.10','NC_000018.10','NC_000020.11','NC_000022.11']
chr_list_length=len(chr_list)

prediction={}
prediction_count={}
for the_chr in chr_list:
    input_data=pickle.load(open('../../data/'+the_chr+'.pkl', "rb"))
    length_chr=input_data.shape[0]-1024
    prediction[the_chr]=np.zeros(length_chr)
    prediction_count[the_chr]=np.zeros(length_chr)+0.0001
    i=0
    i_ccc=i
    chunk_i=0
    all_sample=[]
    while (i<length_chr):
        sample=input_data[i:i+1024,:]
        all_sample.append(sample)
        chunk_i=chunk_i+1
        i=i+256
        if ((chunk_i%2000)==0):
            print(np.asarray(all_sample).shape)
            pred0=np.asarray(model0.predict(np.asarray(all_sample)))
            pred1=np.asarray(model1.predict(np.asarray(all_sample)))
            pred2=np.asarray(model2.predict(np.asarray(all_sample)))
            pred3=np.asarray(model3.predict(np.asarray(all_sample)))
            pred4=np.asarray(model4.predict(np.asarray(all_sample)))
            pred0=pred0.reshape(pred0.shape[0],pred0.shape[1])
            pred1=pred1.reshape(pred1.shape[0],pred1.shape[1])
            pred2=pred2.reshape(pred2.shape[0],pred2.shape[1])
            pred3=pred3.reshape(pred3.shape[0],pred3.shape[1])
            pred4=pred4.reshape(pred4.shape[0],pred4.shape[1])
            pred=pred0+pred1+pred2+pred3+pred4
            pred=pred/5 # n*1024
            iii=0
            while (iii<pred.shape[0]):
                prediction[the_chr][i_ccc+256:i_ccc+768]=pred[iii][256:768]+prediction[the_chr][i_ccc+256:i_ccc+768]
                prediction_count[the_chr][i_ccc+256:i_ccc+768]=prediction[the_chr][i_ccc+256:i_ccc+768]+1
                iii=iii+1
                i_ccc=i_ccc+256
            all_sample=[]
            chunk_i=0
            i_ccc=i
    #        print('done chunk')
    #        break
    prediction[the_chr]=np.divide(prediction[the_chr],prediction_count[the_chr])
    np.save(the_chr+'.pred',prediction[the_chr])
