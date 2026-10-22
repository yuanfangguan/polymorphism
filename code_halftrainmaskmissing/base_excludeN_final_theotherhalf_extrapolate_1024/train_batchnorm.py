import cnn_batchnorm as cnn
import keras
import numpy as np
import tensorflow as tf
import pickle
import random

random.seed(3.14)
model=cnn.cnn1d(1024,4)


#chr_list=['NC_000002.12','NC_000004.12','NC_000006.12','NC_000008.11','NC_000010.11','NC_000012.12','NC_000014.9','NC_000016.10','NC_000018.10','NC_000020.11','NC_000022.11']
chr_list=['NC_000004.12','NC_000006.12','NC_000008.11','NC_000010.11','NC_000012.12','NC_000014.9','NC_000016.10','NC_000018.10','NC_000020.11','NC_000022.11']
chr_list_length=len(chr_list)

all_data={}
all_label={}
total_len_train={}
total_len_val={}
for the_chr in chr_list:
    all_data[the_chr]=pickle.load(open('../../data/'+the_chr+'.pkl', "rb"))
    all_label[the_chr]=pickle.load(open('../../data/'+the_chr+'_label.pkl', "rb"))
    total_len_train[the_chr]=all_data[the_chr].shape[0]-1024


batch_size=16

def generate_data(data_input, label_input, batch_size,if_train):
    while True:
        data_batch=[] # features to be fed in
        label_batch=[] # corresponding gold standard
        for aaa in range(batch_size):
            chr_id=int(random.random()*chr_list_length)
            the_chr=chr_list[chr_id]
            start=int(random.random()*total_len_train[the_chr])
            sample=data_input[the_chr][start:start+1024,:] 
            label=label_input[the_chr][start:start+1024] # minimal 0, maximal 1024
            if (np.min(np.max(sample,axis=1))==0):
                pass
            if (np.sum(label)==0):
                pass
            else:
                data_batch.append(sample)
                label_batch.append(label)
        data_batch=np.asarray(data_batch)
        label_batch=np.asarray(label_batch)
        yield data_batch,label_batch

import sys
callbacks=[
   keras.callbacks.ModelCheckpoint('model'+sys.argv[1],
   verbose=0,save_weights_only=False,save_best_only=True, monitor='val_loss')

]
from tensorflow.keras.callbacks import EarlyStopping
early_stopping = EarlyStopping(monitor='val_loss',patience=2)

model.fit_generator(
    generate_data(all_data, all_label, batch_size, True),
    #steps_per_epoch=int(data_input.shape[2]//batch_size),nb_epoch=10,
    steps_per_epoch=1000,epochs=50, # 500*16*10 examples passing by
    validation_data=generate_data(all_data,all_label, batch_size,False), ## not a strict partition but due to the large size of DNA we can safely use the random draw from the same sample pool.
    validation_steps=20000, # 2000*16 examples= 32000-- we are not to evaluate obviously the whole DNA
    callbacks=callbacks)
