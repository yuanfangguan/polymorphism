import keras
import tensorflow as tf
import numpy as np
from keras.models import Model
from keras.layers import Input, Conv1D, MaxPooling1D, BatchNormalization, Dense, Flatten,Layer,Embedding,Bidirectional,CuDNNLSTM,Dropout,TimeDistributed,Lambda,Activation,GlobalAveragePooling1D, Concatenate,concatenate,Conv1DTranspose
#from keras.optimizers import Adam
from tensorflow.keras.optimizers import Adam
from keras.initializers import glorot_uniform,Ones, Zeros
from keras.losses import binary_crossentropy
from keras import backend as K
from keras import losses



## same as cross entropy
def self_crossentropy(y_true, y_pred):
    y_true = K.flatten(y_true)
    y_pred = tf.clip_by_value(K.flatten(y_pred), 1e-7, (1.0 - 1e-7))
    out = -(y_true * K.log(y_pred) + (1.0 - y_true) * K.log(1.0 - y_pred))
    return K.mean(out)


# mask out the part you don't want
def crossentropy_cut(y_true,y_pred):
    y_true_f = K.flatten(y_true)
    y_pred_f = K.flatten(y_pred)
    y_pred_f= tf.clip_by_value(y_pred_f, 1e-7, (1. - 1e-7))
    mask=K.cast(K.greater_equal(y_true_f,-0.5),dtype='float32')
    out = -(y_true_f * K.log(y_pred_f)*mask + (1.0 - y_true_f) * K.log(1.0 - y_pred_f)*mask)
    out=K.mean(out)
    return out


# same as mean squared error
def mean_squared_error(y_true, y_pred):
    return K.mean(K.square(y_pred - y_true), axis=-1)

# clip at a range
def clipped_mse(y_true, y_pred):
    return K.mean(K.square(K.clip(y_pred, 0., 1.0) - K.clip(y_true, 0., 1.0)), axis=-1)



def cnn1d(size, channel): ## size example 2048*4

    the_input= Input((size,channel)) #2048*4

    conv1= BatchNormalization()(Conv1D(16, kernel_size = 3, activation='relu',padding='same')(the_input))#2048 * 16
    conv1= BatchNormalization()(Conv1D(16, kernel_size = 3, activation='relu',padding='same')(conv1))
    maxpool1= MaxPooling1D(pool_size=2)(conv1) #1024*16
    conv2= BatchNormalization()(Conv1D(32, kernel_size = 3, activation='relu',padding='same')(maxpool1)) #1024*32
    conv2= BatchNormalization()(Conv1D(32, kernel_size = 3, activation='relu',padding='same')(conv2))
    maxpool2= MaxPooling1D(pool_size=2)(conv2) # 512*32
    conv3= BatchNormalization()(Conv1D(64, kernel_size = 3, activation='relu',padding='same')(maxpool2)) #512*64
    conv3= BatchNormalization()(Conv1D(64, kernel_size = 3, activation='relu',padding='same')(conv3))
    maxpool3= MaxPooling1D(pool_size=2)(conv3) # 256*64
    conv4= BatchNormalization()(Conv1D(128, kernel_size = 3, activation='relu',padding='same')(maxpool3)) #256*128
    conv4= BatchNormalization()(Conv1D(128, kernel_size = 3, activation='relu',padding='same')(conv4))
    maxpool4= MaxPooling1D(pool_size=2)(conv4) # 128*128 
    conv5= BatchNormalization()(Conv1D(256, kernel_size = 3, activation='relu',padding='same')(maxpool4)) #128*256
    conv5= BatchNormalization()(Conv1D(256, kernel_size = 3, activation='relu',padding='same')(conv5)) #128*256

    up6 = concatenate([Conv1DTranspose(256, 2, strides=2, padding='same')(conv5), conv4], axis=2) # conv3 size
    conv6= Conv1D(256, 3, activation='relu', padding='same')(up6)
    conv6= Conv1D(256, 3, activation='relu', padding='same')(conv6)
    up7 = concatenate([Conv1DTranspose(128, 2, strides=2, padding='same')(conv6), conv3], axis=2) # conv2 size
    conv7= Conv1D(128, 3, activation='relu', padding='same')(up7)
    conv7= Conv1D(128, 3, activation='relu', padding='same')(conv7)
    up8 = concatenate([Conv1DTranspose(64, 2, strides=2, padding='same')(conv7), conv2], axis=2) # conv2 size
    conv8= Conv1D(64, 3, activation='relu', padding='same')(up8)
    conv8= Conv1D(64, 3, activation='relu', padding='same')(conv8)
    up9 = concatenate([Conv1DTranspose(32, 2, strides=2, padding='same')(conv8), conv1], axis=2) # conv2 size
    conv9= Conv1D(32, 3, activation='relu', padding='same')(up9)
    conv9= Conv1D(1, 3, activation='sigmoid', padding='same')(conv9)

    #y= Dense(1)(layer) # 1 refers to how many things you want to predict 

    model = Model(inputs=[the_input], outputs=[conv9])
    model.compile(loss=self_crossentropy, optimizer='adam',metrics = ['mean_absolute_error'])
    #model.compile(loss=mean_squared_error, optimizer='adam',metrics = ['mean_absolute_error'])
   # model.compile(loss='BinaryCrossentropy', lr=0.0001, optimizer='adam',metrics = ['mean_absolute_error'])

    return model
