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
import matplotlib

matplotlib.rcParams['mathtext.fontset'] = 'custom'
matplotlib.rcParams['mathtext.rm'] = 'Bitstream Vera Sans'
matplotlib.rcParams['mathtext.it'] = 'Bitstream Vera Sans:italic'
matplotlib.rcParams['mathtext.bf'] = 'Bitstream Vera Sans:bold'
#matplotlib.pyplot.title(r'ABC123 vs $\mathrm{ABC123}^{123}$')

matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
matplotlib.rcParams.update({'font.size': 16})
#matplotlib.pyplot.title(r'ABC123 vs $\mathrm{ABC123}^{123}$')

file = open('FinalDF_50.pckl', 'rb')
DF = pickle.load(file)
DF = pd.DataFrame(DF)
DF = DF.dropna(how='all')
print(DF.shape)
print(type(DF))

DF_e1 = DF[DF['Early Finish2'] == False]
print(DF_e1.shape)
print((len(DF.index)- len(DF_e1.index))/ len(DF.index))

decomp = list(DF['Heat of Decomposition'])
new_decomp = list(DF_e1['Heat of Decomposition'])

for i, v in enumerate(new_decomp):
    if v > 0:
        new_decomp[i] = 0

for i, v in enumerate(decomp):
    if v > 0:
        decomp[i] = 0


print(sum(new_decomp)/len(new_decomp))
print(sum(DF_e1['e_above_hull'])/len(DF_e1['e_above_hull']))

abs_difference = abs(-np.array(decomp) - np.array(DF['e_above_hull']))
print("with early finshers (MAE): ", sum(abs_difference)/len(abs_difference))

new_difference = (-np.array(new_decomp) - np.array(DF_e1['e_above_hull']))
print(len([d for i, d  in enumerate(new_difference) if d>0.001]))
print([DF.iloc[i] for i, d  in enumerate(new_difference) if d>0.001 and DF['Formula'][i] == 'CrCl3'])

print(np.array(new_decomp), np.array(DF_e1['e_above_hull']))
print("without early finshers: ", sum(new_difference)/len(new_difference))

sqr_difference = (-np.array(decomp) - np.array(DF['e_above_hull']))**2
print("RMSE: ", math.sqrt(sum(sqr_difference)/len(sqr_difference)))

'''
# MEAN ABSOLUTE ERROR
fig, axes = plt.subplots()
axes.hist(abs_difference, bins=6000, color = 'blue', edgecolor = 'black')
axes.set_xlim(-0.01, 0.05)
axes.set_xlabel(r'$|E_{greedy} - E_{hull}|$ [eV/atom]', fontsize=18)
axes.set_ylabel('Frequency',fontsize=16)
#plt.tight_layout()
plt.gcf().subplots_adjust(bottom=0.15)
plt.style.use('classic')
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0), useMathText=True)
#axes.tick_params(axis="y", direction='in', which='both')
#axes.tick_params(axis="x", direction='in', which='both')
#axes.tick_params(bottom=True, top=True, left=True, right=True)
plt.show()
#plt.savefig('../Downloads/fig.pdf')
'''

print(len([i for i in new_decomp if i == 0]))
print(len([i for i in DF['e_above_hull'] if i ==0]))
print(max(DF['e_above_hull']))

'''
zero = 0
non_zero = 0
for i, v in enumerate(new_decomp):
    if v == 0 and hull[i] == 0
'''
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
'''
for i, l in enumerate(DF['elements']):
    if l == ['Ti']:
        print(DF.iloc[i])
'''
