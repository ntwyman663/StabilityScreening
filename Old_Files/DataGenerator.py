#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  3 19:10:25 2019

@author: NickT
"""

import pickle
from pymatgen import MPRester

mat_api_key = 'JWV6Fi4f6VfxROtHO2uP'

mpr = MPRester(mat_api_key)

print("Loading Compounds....")
Querey = mpr.query(criteria = {'elements': ['Si'], 'nelements': 1}, properties=["task_id", "pretty_formula", 'e_above_hull',
                          'elements', 'volume', 'formation_energy_per_atom', 'band_gap',
                           'nsites', 'unit_cell_formula'])

file = open('MPDatabase.pickle', 'wb')
pickle.dump(all_compounds, file)
file.close()



stable_phase = []

for compound in all_compounds: #find all compounds with e_above_hull within 0.05 of 0
    if abs(compound['e_above_hull']) < criteria/1000:
        stable_phase.append(compound)

