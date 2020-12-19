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

file = open('FinalDF_50.pckl', 'rb')
DF = pickle.load(file)
DF = DF.dropna(how='all')
print(DF.shape)
#
#file = open('FinalDF_50.0.1.pckl', 'rb')
#DF = pickle.load(file)
#
#file2 = open('num_elements50.0.pckl', 'rb')
#num_elements = pickle.load(file2)


#x = [random.gauss(40,20) for _ in range(400)]
#y = [random.gauss(4,2) for _ in range(400)]
#d = {'col1': x, 'col2': y}
#df = pd.DataFrame(data=d)
#
#sns.distplot(df['col1'], color = 'g')
#plt.xlabel('categories')
#plt.ylabel('values')
#
#sns.plt.show()

# DIFFERENCE BETWEEN HEAT OF DECOMP AND HULL
difference = list((-DF['Heat of Decomposition'] - DF['e_above_hull']))
sns.distplot(difference, color = 'g')
plt.xlabel('Heat of Decomposition Error')
plt.ylabel('Normalised Frequency')
plt.show()

plt.boxplot(difference)
plt.show()
print(sum(difference)/len(difference))





'''
#NUMBER OF COMPETING PHASES
x = [i for i in DF['Competing Phase Number (with formation E correction)'] if not math.isnan(i) ]
sns.distplot(x, bins = 500, color = 'g')
plt.xlim(0, 1200)
plt.xlabel('Number of Competing Phases')
plt.ylabel('Normalised Frequency')
plt.show()

#NUMBER OF COMPETING PHASES 0-100
print(len(DF['Competing Phase Number (with formation E correction)']))
x = [i for i in DF['Competing Phase Number (with formation E correction)'] if not math.isnan(i) and i < 101 ]
#print(len(nans))
sns.distplot(x, color = 'g')
plt.xlim(0, 100)
plt.xlabel('Number of Competing Phases')
plt.ylabel('Normalised Frequency')
plt.show()

#NUMBER OF ELEMENTS IN 
norm = 0
for i in num_elements: 
    if i ==3: 
        norm +=1

sns.countplot(num_elements, color = 'g')
#plt.xlim(0, 100)
plt.xlabel('Number of Elements in a Material')
plt.ylabel('Normalised Frequency')
plt.xlim(0, 5)
plt.rc('xtick', labelsize=16) 
plt.rc('ytick', labelsize=16) 
plt.show()


#NUM ELEMENTS IN VS. NUM COMPETING PHASES
x = [i for i in DF['Competing Phase Number (with formation E correction)'] if not math.isnan(i)]
y = [x[i] for i in range(len(x)) if x[i]<31]
z = [num_elements[i] for i in range(len(x)) if x[i]<31]

y2 = [y[i] for i in range(len(y)) if z[i]<7]
z2 = [z[i] for i in range(len(y)) if z[i]<7]
DList = {'Number of Competing Phases':y2, 'Number of Elements in a Material': z2}
df1 = pd.DataFrame.from_dict(DList)
sns.jointplot(x="Number of Competing Phases", y="Number of Elements in a Material", 
              data=df1, kind="kde", color = 'g')


'''




#NUMBER OF POTENTIAL
x = [i for i in DF['Number of Oxides (with formation E correction)'] if not math.isnan(i) ]
sns.distplot(x, bins = 500, color = 'r')
plt.xlim(0, 1200)
plt.xlabel('Number of Potential Oxides')
plt.ylabel('Normalised Frequency')
plt.show()


#NUMBER OF POTENTIAL OXIDES 0-400
x = [i for i in DF['Number of Oxides (with formation E correction)'] if not math.isnan(i) and i < 401 ]
sns.distplot(x, color = 'r')
plt.xlim(0, 400)
plt.xlabel('Number of Potential Oxides')
plt.ylabel('Normalised Frequency')
plt.show()

'''
#NUMBER OF ELEMENTS IN 

sns.distplot(num_elements, color = 'r')
#plt.xlim(0, 100)
plt.xlabel('Number of Elements in a Material')
plt.ylabel('Normalised Frequency')
plt.show()


#NUM ELEMENTS IN VS. NUM COMPETING PHASES
x = [i for i in DF['Number of Oxides (with formation E correction)'] if not math.isnan(i)]
y = [x[i] for i in range(len(x)) if x[i]<301]
z = [num_elements[i] for i in range(len(x)) if x[i]<301]

y2 = [y[i] for i in range(len(y)) if z[i]<7]
z2 = [z[i] for i in range(len(y)) if z[i]<7]
DList = {'Number of Potential Oxides':y2, 'Number of Elements in a Material': z2}
df1 = pd.DataFrame.from_dict(DList)
sns.jointplot(x="Number of Potential Oxides", y="Number of Elements in a Material", 
              data=df1, kind="kde", color = 'r')

'''
#HEAT OF DECOMPOSITION 
x = [i for i in DF['Heat of Decomposition'] if not math.isnan(i) and i>-0.1]
zeros = [i for i in x if i == 0]
#print(len(x), len(zeros))
sns.distplot(x, color = 'g', bins=500)
plt.xlim(-0.1, 0)
plt.xlabel('Heat of Decomposition /eV/atom')
plt.ylabel('Normalised Frequency')
plt.rc('xtick', labelsize=12) 
plt.rc('ytick', labelsize=12) 
plt.show()


"""
#HEAT OF DECOMPOSITION WITHOUT EARLY FINISHERS
z = [i for i in DF['Heat of Decomposition'] if not math.isnan(i)]
y = [i for i in DF['Early Finish1'] if not math.isnan(i)]
x2 = [z[i] for i in range(len(y)) if y[i]==False]
sns.distplot(x2, color = 'g')
plt.xlim(-0.1, 0)
plt.xlabel('Heat of Decomposition /eV/atom')
plt.ylabel('Normalised Frequency')
plt.rc('xtick', labelsize=12) 
plt.rc('ytick', labelsize=12) 
plt.show()
"""

#COMPLEMENTARY GROUP SIZE OF COMPLETED
z = [i for i in DF['Number of Complementary Phases'] if not math.isnan(i)]
y = [i for i in DF['Early Finish1'] if not math.isnan(i)]
x3 = [z[i] for i in range(len(y)) if y[i]==False]
sns.distplot(x3, color = 'g')
plt.xlabel('Heat of Decomposition /eV/atom')
plt.ylabel('Normalised Frequency')
plt.xlim(0, 7)
plt.rc('xtick', labelsize=12) 
plt.rc('ytick', labelsize=12) 
plt.show()

'''
x4 = [num_elements[i] for i in range(len(y)) if y[i]==False]

DList = {'Heat of Decomposition /eV/atom':x2, 'Number of Complementary Phases': x4}
df1 = pd.DataFrame.from_dict(DList)
sns.jointplot(x='Heat of Decomposition /eV/atom', y='Number of Complementary Phases', 
              data=df1, kind="kde", color = 'g')
'''



##HEAT OF OXIDATION IGNORING ZERO
x = [i for i in DF['Complementary Heat of Oxidation'] if not math.isnan(i)]
sns.distplot(x, color = 'r')
#plt.xlabel('Number of Competing Phases')
#plt.ylabel('Normalised Frequency')
#plt.xlim(-12, 0)
#sns.plt.show()

#HEAT OF OXIDATION WITHOUT EARLY FINISHERS
z = [i for i in DF['Heat of Oxidation'] if not math.isnan(i)]
y = [i for i in DF['Early Finish2'] if not math.isnan(i)]
x2 = [z[i] for i in range(len(y)) if y[i]==False]
sns.distplot(x2, color = 'c')
plt.xlim(-12, 0)
plt.xlabel('Heat of Oxidation /eV/atom')
plt.ylabel('Normalised Frequency')
plt.rc('xtick', labelsize=12) 
plt.rc('ytick', labelsize=12)
plt.show()


#HEAT OF OXIDATION WITHOUT EARLY FINISHERS and 0.1 constraint
z = [i for i in DF['Heat of Oxidation'] if not math.isnan(i)]
y = [i for i in DF['Early Finish2'] if not math.isnan(i)]
x2 = [z[i] for i in range(len(y)) if y[i]==False and z[i] <= -0.1]
sns.distplot(x2, color = 'c')
plt.xlim(-12, 0)
plt.xlabel('Heat of Oxidation /eV/atom')
plt.ylabel('Normalised Frequency')
plt.rc('xtick', labelsize=12) 
plt.rc('ytick', labelsize=12)
plt.show()



#HEAT OF OXIDATION WITHOUT EARLY FINISHERS 2
z = [i for i in DF['Complementary Heat of Oxidation'] if not math.isnan(i)]
y = [i for i in DF['Early Finish2'] if not math.isnan(i)]
x2 = [z[i] for i in range(len(y)) if y[i]==False]
sns.distplot(x2, color = 'c')
plt.xlim(-12, 0)
plt.xlabel('Heat of Decomposition /eV/atom')
plt.ylabel('Normalised Frequency')
plt.rc('xtick', labelsize=12) 
plt.rc('ytick', labelsize=12) 
plt.show()

print((np.mean(x2), np.std(x2)))

#COMPLEMENTARY GROUP SIZE OF COMPLETED
z = [i for i in DF['Number of Complementary Oxides'] if not math.isnan(i)]
y = [i for i in DF['Early Finish1'] if not math.isnan(i)]
x2 = [z[i] for i in range(len(y)) if y[i]==False]
sns.countplot(x2, color = 'g')
plt.xlabel('Heat of Decomposition /eV/atom')
plt.ylabel('Normalised Frequency')
plt.xlim(0, 7)
plt.rc('xtick', labelsize=12) 
plt.rc('ytick', labelsize=12) 
plt.show()


#BEST VOLUME RATIO
x = [i for i in DF['ID of Best Volume Ratio'] if not math.isnan(i) and i < 3]
sns.distplot(x, color = 'b')

plt.xlabel('Optimal PBR')
plt.ylabel('Normalised Frequency')
plt.rc('xtick', labelsize=12) 
plt.rc('ytick', labelsize=12) 
plt.xlim(0, 3)
plt.show()

z = [i for i in range(len(DF['Formula'])) if DF['Formula'][i] == 'Cu']
m = [DF['ID of Best Volume Ratio'][i] for i in z]
m2 = [DF['Best Volume Ratio'][i] for i in z]

r = [i for i in range(len(DF['ID of Best Volume Ratio'])) if not math.isnan(DF['ID of Best Volume Ratio'][i]) and -0.35 <=1-DF['ID of Best Volume Ratio'][i]<=0]
r1 = [DF['ID of Best Volume Ratio'][i] for i in r]
r2 = [DF['Formula'][i] for i in r]
r3 = DF['ID of Best Volume Ratio'][:]

#plt.legend()
#y = [i for i in DF['Early Finish2'] if not math.isnan(i)]
#z = [x[i] for i in range(len(y)) if y[i]==True]
