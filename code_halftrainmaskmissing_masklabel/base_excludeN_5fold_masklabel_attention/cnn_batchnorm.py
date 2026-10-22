import keras
from tensorflow.keras import layers, models

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

def attention_gate(skip_connection, gating_signal, inter_channels):
    """
    Implements a simple attention gate for 1D input.
    """
    theta_x = layers.Conv1D(inter_channels, kernel_size=1, strides=1, padding='same')(skip_connection)
    phi_g = layers.Conv1D(inter_channels, kernel_size=1, strides=1, padding='same')(gating_signal)
    add = layers.Add()([theta_x, phi_g])
    relu = layers.ReLU()(add)
    psi = layers.Conv1D(1, kernel_size=1, strides=1, padding='same')(relu)
    sigmoid = layers.Activation('sigmoid')(psi)
    return layers.Multiply()([skip_connection, sigmoid])

def unet_with_attention(input_length):
    """
    Defines a U-Net-like architecture with attention gates for 1D sequences.
    """
    inputs = layers.Input(shape=(input_length, 4))  # Assuming input has 4 features per position

    # Encoder
    c1 = layers.Conv1D(64, kernel_size=3, activation='relu', padding='same')(inputs)
    p1 = layers.MaxPooling1D(pool_size=2)(c1)

    c2 = layers.Conv1D(128, kernel_size=3, activation='relu', padding='same')(p1)
    p2 = layers.MaxPooling1D(pool_size=2)(c2)

    c3 = layers.Conv1D(256, kernel_size=3, activation='relu', padding='same')(p2)
    p3 = layers.MaxPooling1D(pool_size=2)(c3)

    # Bottleneck
    c4 = layers.Conv1D(512, kernel_size=3, activation='relu', padding='same')(p3)

    # Decoder with Attention Gates
    u3 = layers.UpSampling1D(size=2)(c4)
    a3 = attention_gate(c3, u3, inter_channels=256)
    d3 = layers.Conv1D(256, kernel_size=3, activation='relu', padding='same')(u3)
    d3 = layers.Concatenate()([d3, a3])

    u2 = layers.UpSampling1D(size=2)(d3)
    a2 = attention_gate(c2, u2, inter_channels=128)
    d2 = layers.Conv1D(128, kernel_size=3, activation='relu', padding='same')(u2)
    d2 = layers.Concatenate()([d2, a2])

    u1 = layers.UpSampling1D(size=2)(d2)
    a1 = attention_gate(c1, u1, inter_channels=64)
    d1 = layers.Conv1D(64, kernel_size=3, activation='relu', padding='same')(u1)
    d1 = layers.Concatenate()([d1, a1])

    # Output Layer
    outputs = layers.Conv1D(1, kernel_size=1, activation='sigmoid')(d1)
    model = Model(inputs=[inputs], outputs=[outputs])
    model.compile(loss=self_crossentropy, optimizer='adam',metrics = ['mean_absolute_error'])

    return model

