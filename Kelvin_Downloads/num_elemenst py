#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr  7 09:58:28 2019

@author: NickT
"""

from pymatgen import MPRester
import pandas as pd
import matplotlib.pyplot as plt
import csv
import multiprocessing as mp
import pickle
import tqdm
import time

mat_api_key = 'JWV6Fi4f6VfxROtHO2uP'

mpr = MPRester(mat_api_key)

print("Loading Compounds....")
file = open('MPDatabase.pickle', 'rb')
all_compounds = pickle.load(file)

criteria = float(input("Enter Stable Phase Criteria in meV: ")) 

print('Finding Stable Phases....')    
    



num_elements = []
for compound in all_compounds:
    if abs(compound['e_above_hull']) < criteria/1000: #if stable
        num_elements.append(len(compound['elements']))


filename = 'num_elements' + str(criteria) + '.pckl'
f = open(filename, 'wb')
pickle.dump(num_elements, f)
f.close()
