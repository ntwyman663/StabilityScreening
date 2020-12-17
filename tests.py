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


DF_e1 = DF[DF['Early Finish2'] == True]
print(DF_e1.shape)
