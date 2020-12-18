#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  4 20:12:37 2019

@author: NickT
"""
import pickle
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import random
import seaborn as sns
from scipy import stats
import math

file = open('decomposition_test.pckl', 'rb')
DF = pickle.load(file)
DF = pd.DataFrame(DF)
DF = DF.dropna(how='all')
print(DF.shape)


#DF_e1 = DF[DF['Early Finish2'] == True]
#print(DF_e1.shape)

new_decomp = list(DF['Heat of Decomposition'])
for i, v in enumerate(new_decomp):
    if v > 0:
        new_decomp[i] = 0

print(sum(new_decomp)/len(new_decomp))
print(sum(DF['e_above_hull'])/len(DF['e_above_hull']))

difference = (-pd.Series(new_decomp) - DF['e_above_hull'])
print(sum(difference)/len(difference))

plt.boxplot(difference)
plt.show()

print(len([i for i in new_decomp if i == 0]))
print(len([i for i in DF['e_above_hull'] if i ==0]))
print(max(DF['e_above_hull']))

'''
zero = 0
non_zero = 0
for i, v in enumerate(new_decomp):
    if v == 0 and hull[i] == 0
'''
