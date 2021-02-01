# -*- coding: utf-8 -*-
"""
Created on Sun Jan 31 16:27:37 2021

@author: Meenal
"""

import functools
import itertools
import os
import random

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from sklearn import metrics
from sklearn import preprocessing
from sklearn.model_selection import train_test_split, KFold, cross_val_score
from keras.callbacks import ModelCheckpoint
from keras.layers import Dense, Activation, Flatten, Dropout, Reshape
from keras.layers import Conv1D,Conv2D, MaxPooling2D
from keras.models import Sequential,Model
from keras.utils.np_utils import to_categorical
from keras import optimizers
from keras.optimizers import Adam,SGD
from keras.layers.normalization import BatchNormalization
from keras.regularizers import l2
import copy
from sklearn.metrics import roc_auc_score
from Bio import SeqIO
import pandas as pd
from keras.layers import Bidirectional,TimeDistributed
import keras
from keras.models import Sequential
from keras.layers import Dense, Dropout, Flatten, Reshape, Lambda, LSTM
from keras.layers import Conv2D, MaxPooling2D
from keras.callbacks import ModelCheckpoint, ReduceLROnPlateau
import numpy as np
from Bio import SeqIO
from numpy import array
from numpy import argmax
from sklearn.model_selection import train_test_split
from imblearn.under_sampling import RandomUnderSampler
from sklearn.utils import shuffle
from keras.layers.embeddings import Embedding
from keras import backend as K
from keras.backend import expand_dims
import matplotlib.pyplot as plt
from keras.regularizers import l1, l2
from sklearn.metrics import roc_curve, auc, classification_report
from keras.models import load_model
import pandas as pd
import time
res = 'ST' #'Y'
r_test_x = []
r_test_y = []

# Test DATASET -------------------------------------------------------------
#for positive sequence
def inner1():
    #Input
    data = seq_record.seq
    #data = data[cut_off:-cut_off]
    #rint(data) 
    # integer encode input data
    for char in data:
        if char not in alphabet:
            return
    integer_encoded = [char_to_int[char] for char in data]
    r_test_x.append(integer_encoded)
    r_test_y.append(posit_1)
for seq_record in SeqIO.parse("dataset/test_Pos_ST.fasta", "fasta"):
    inner1()
# for negative sequence
def inner2():
    #Input
    data = seq_record.seq
    #data = data[cut_off:-cut_off]
    #print(data) 
    # integer encode input data
    for char in data:
        if char not in alphabet:
            return
    integer_encoded = [char_to_int[char] for char in data]
    r_test_x.append(integer_encoded)
    r_test_y.append(negat_0)
for seq_record in SeqIO.parse("dataset/test_Neg_ST.fasta", "fasta"):
    inner2()
# Changing to array (matrix)    
r_test_x = array(r_test_x)
r_test_y1 = array(r_test_y)#pd.read_csv("traintest/OldDephos_labels_train_S_33_clean.csv",header= [0])#array(train_y)
r_test_y = keras.utils.to_categorical(r_test_y1, num_classes)

model = load_model('ComDephos_ST.h5')#+str(res)+'
score = model.evaluate(x_test,y_test, verbose=0)
print('Train-val loss:', score[0])
print('Train-val accuracy:', score[1])
acc_train = score[1]
from sklearn.metrics import matthews_corrcoef
from sklearn.metrics import confusion_matrix
Y_pred = model.predict(x_test)
t_pred2 = Y_pred[:,1]
Y_pred = (Y_pred > 0.5)
y_pred1 = [np.argmax(y, axis=None, out=None) for y in Y_pred]
y_pred1 = np.array(y_pred1)

print("Matthews Correlation : ",matthews_corrcoef(y_test[:,1],y_pred1))
print("Confusion Matrix : \n",confusion_matrix( y_test[:,1],y_pred1))
mcc_train = matthews_corrcoef(y_test[:,1],y_pred1)
# For sensitivity and specificity
sp_1, sn_1 = confusion_matrix(y_test[:,1], y_pred1)
sp_2_train = sp_1[0]/(sp_1[0]+sp_1[1])
sn_2_train = sn_1[1]/(sn_1[0]+sn_1[1])
# ROC

fpr, tpr, _ = roc_curve(y_test[:,1], t_pred2)
roc_auc_train = auc(fpr, tpr)
print("AUC : ", roc_auc_train)
print(classification_report(y_test[:,1], y_pred1))
print("Specificity = ",sp_2_train, " Sensitivity = ",sn_2_train)
plt.figure()
lw = 2
plt.plot(fpr, tpr, color='darkorange',
         lw=lw, label='ROC curve (area = %0.2f)' % roc_auc_train)
plt.plot([0, 1], [0, 1], color='navy', lw=lw, linestyle='--')
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('ROC curve for'+str(ty))
plt.legend(loc="lower right")
plt.show()


    