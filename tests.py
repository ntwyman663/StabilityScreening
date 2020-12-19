#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: NickT

This file contains tests to analyse teh .pckl file produced by Analyser.py
"""
import pickle
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import random
import seaborn as sns
from scipy import stats
import math

file = open('FinalDF_50.pckl', 'rb')
DF = pickle.load(file)
DF = pd.DataFrame(DF)
DF = DF.dropna(how='all')
print(DF.shape)
print(type(DF))

#DF_e1 = DF[DF['Early Finish2'] == True]
#print(DF_e1.shape)
'''
new_decomp = list(DF['Heat of Decomposition'])
for i, v in enumerate(new_decomp):
    if v > 0:
        new_decomp[i] = 0

print(sum(new_decomp)/len(new_decomp))
print(sum(DF['e_above_hull'])/len(DF['e_above_hull']))

difference = (-pd.Series(new_decomp) - DF['e_above_hull'])
print(sum(difference)/len(difference))

plt.boxplot(difference)
plt.ylim(-0.2, 0.2)
plt.show()
sns.distplot(difference, bins=500)
plt.xlim(-0.2, 0.2)
plt.show()

print(len([i for i in new_decomp if i == 0]))
print(len([i for i in DF['e_above_hull'] if i ==0]))
print(max(DF['e_above_hull']))
'''
'''
zero = 0
non_zero = 0
for i, v in enumerate(new_decomp):
    if v == 0 and hull[i] == 0
'''

def pbr(val):
    return val>1 and val<2

first = []
for i, l in enumerate(DF['Volume Ratios']):
    if len(l) > 1:
        if not DF['Early Finish2'][i] and pbr(l[1]) and pbr(l[0]) and \
                DF['Heat of Oxidation'][i] < 0:
            first.append(DF.iloc[i])
        #elif not DF['Early Finish2'][i] and pbr(l[0]):
            #first.append(DF.iloc[i])
    else:
        if not DF['Early Finish2'][i] and pbr(l[0]) and DF['Heat of Oxidation'][i] < 0:
            first.append(DF.iloc[i])

print(len(first))
singles = []
for i in first:
    if len(i['elements']) == 1: 
        singles.append((i['Formula'], i['Heat of Oxidation']))
#singles = list(dict.fromkeys(singles))
print(sorted(singles, key= lambda x:x[0]))
print(sorted(singles, key= lambda x:x[1]))
'''
for i, l in enumerate(DF['elements']):
    if l == ['Ti']:
        print(DF.iloc[i])
'''
