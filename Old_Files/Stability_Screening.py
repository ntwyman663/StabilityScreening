#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  4 18:07:03 2019

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
all_compounds = mpr.query({}, properties=["task_id", "pretty_formula", 'e_above_hull',
                          'elements', 'volume', 'formation_energy_per_atom', 'band_gap',
                          'unit_cell_formula', 'nsites'])

criteria = float(input("Enter Stable Phase Criteria in meV: ")) 


print('Finding Stable Phases....')    
    
stable_phase = []

for compound in tqdm.tqdm(all_compounds): #find all compounds with e_above_hull within 0.05 of 0
    if abs(compound['e_above_hull']) < criteria/1000:
        stable_phase.append(compound)

#pool = mp.Pool(processes=1)
#
#stable_phase = list(tqdm.tqdm(pool.imap(find_stable_phases, all_compounds), total=86680))


######## COMPETING PHASE AND OXIDE CALCULATION ########

def find_comp(stable_oxides, compound_unit_cell, compound_formE, condition):
    '''
    Finds complementary oxide or competing phases group and associated total heat of oxidation/
    decomposition
    
    args:
        stable_oxides = list of dictionaries of stable oxides or competing phases
        compound_unit_cell = dict of elements in unit cell of original compound
        compound_formE = formation energy of original compound
        condition = string dictating whether it is for comp oxide or comp competing phases
    
    output:
        tuple: (list of dicitionatries of predicted materials, 
        combined formation energy of these materials (with appropriate ratios),
        heat of oxidation/ decomposition, number of materials in complementary group,
        boolean describing where or not the algorithm finished early)
        
    '''
    result = []
    FinishEarly = False
    #what if positive formE
    if condition == 'Oxide':
        sort_oxides = sorted(stable_oxides, key = lambda oxide: (oxide['formation_energy_per_atom']*oxide['nsites'])/
                             (oxide['nsites']-oxide['unit_cell_formula']['O']))
        
    else: sort_oxides = sorted(stable_oxides, key = lambda oxide: (oxide['formation_energy_per_atom']))
    
    sort_oxides1 = sort_oxides[:]
    orig_natoms = sum(compound_unit_cell.values())
    compound_unit_cell1 = dict(compound_unit_cell)
    for i in compound_unit_cell1:
        compound_unit_cell1[i] = compound_unit_cell1[i]/orig_natoms
    total_formE = 0
    while sum(compound_unit_cell1.values()) != 0 and sort_oxides1 != []: #if all atoms in unit cell not yet accounted for
        oxide = sort_oxides1[0]
        
        intersection = list(set(oxide['elements']).intersection(compound_unit_cell1.keys()))
        intersect_rank = {}
        for element in intersection:
            intersect_rank[element] = compound_unit_cell1[element]/(oxide['unit_cell_formula'][element]
            /oxide['nsites'])
            
        limiting_element = min(intersect_rank, key=intersect_rank.get) #find limiting element
        ratio = intersect_rank[limiting_element] 
        used_up_elements = []
        for element in intersection:
            compound_unit_cell1[element] = compound_unit_cell1[element] - (ratio*
            (oxide['unit_cell_formula'][element]/oxide['nsites']))
            if compound_unit_cell1[element] == 0:
                used_up_elements.append(element)
                
        result.append(oxide)
        sort_oxides1.remove(oxide)
        total_formE += oxide['formation_energy_per_atom']*ratio
        
        sort_oxides1 = [oxide for oxide in sort_oxides1 if 
                        len(set(oxide['elements']).intersection(used_up_elements)) == 0]
                        #remove oxides in list which contain used up elements
        
    
    if sort_oxides1 == [] and sum(compound_unit_cell1.values()) != 0:
        FinishEarly = True
                
    return (result, total_formE, total_formE-compound_formE, len(result), FinishEarly)

    #### FOR TESTING FIND_OXIDES
        
#ABCO4 = {'elements': ['A', 'B', 'C', 'O'], 'formation_energy_per_atom': -750, 'nsites':7,
#        'unit_cell_formula':{'A':1, 'B':1, 'C':1, 'O':4}}
#AO = {'elements': ['A', 'O'], 'formation_energy_per_atom': -100, 'nsites':8,
#        'unit_cell_formula':{'A':4, 'O':4}}
#BO2 = {'elements': ['B', 'O'], 'formation_energy_per_atom': -100, 'nsites':6,
#        'unit_cell_formula':{'B':2, 'O':4}}
#C2O = {'elements': ['C', 'O'], 'formation_energy_per_atom': -300, 'nsites':24,
#        'unit_cell_formula':{'C':16, 'O':8}}
#A2BO6 = {'elements': ['A', 'B', 'O'], 'formation_energy_per_atom': -380, 'nsites':9,
#        'unit_cell_formula':{'A':2, 'B':1, 'O':6}}
#A2CO4 = {'elements': ['A', 'B', 'O'], 'formation_energy_per_atom': -620, 'nsites':63,
#        'unit_cell_formula':{'A':18, 'C':9, 'O':36}}
#
#original = {'A':4, 'B':8, 'C':12}
#
#listt = [ABCO4, AO, BO2, C2O, A2BO6, A2CO4]
#                
#find_comp(listt, original, -400, 'Oxide')    
                    
    ####  
      

def Make_Property_Dict(compound):
    '''
    Function to be iterated over all compounds.
    '''
    PDict = {}
    
    global stable_phase
    
    if abs(compound['e_above_hull']) < criteria/1000: #if stable 
        
        #### FOR NUM PHASES
        competing_phases_id_withform1 = []
        competing_phase_no1 = 0
        comp_listdict =[]
        
        
        #### FOR NUM OXIDES
        v_ratio2 = 0
        oxide_no1 = 0
        oxides_id_withform1 = []
        v_ratio_id2 = 'n/a'
        oxide_listdict = []
            
        elements = compound['elements']
        
        for i in stable_phase:
            #### FOR NUM PHASES
            if set(i['elements']).issubset(elements):
                comp_listdict.append(i) #for find_comp
                
                if i['formation_energy_per_atom'] < compound['formation_energy_per_atom']:
                    #find all other phases containing just those elements                
                    competing_phase_no1 +=1
                    competing_phases_id_withform1.append(i['task_id'])
                    
                    
            #### FOR NUM OXIDES
            if 'O' in i['elements']:
                O = i['unit_cell_formula']['O']
                el = i['elements'][:]
                el.remove('O')
                if set(el).issubset(elements) and len(i['unit_cell_formula']) > 1:
                    oxide_listdict.append(i) #for find_comp
                    
                    if i['formation_energy_per_atom'] < compound['formation_energy_per_atom']:
                        oxide_no1 += 1
                        oxides_id_withform1.append(i['task_id'])
                        
        
        
        
        #### FOR NUM PHASES
        
        PDict['task_id'] = compound['task_id']
        PDict['Formula'] = compound['pretty_formula']
        PDict['Bandgap /eV'] = compound['band_gap']
        
        PDict['Competing Phase Number (with formation E correction)'] = competing_phase_no1
        PDict['Competing Phase List (with formation E correction)'] = competing_phases_id_withform1
        
        y = find_comp(comp_listdict, compound['unit_cell_formula'], compound['formation_energy_per_atom'], 'NotOx')
        PDict['Complementary Competing Phase List'] = y[0]
        PDict['Complementary Heat of Decomposition'] = y[1]
        PDict['Lower Formation Energy Than Original Material'] = y[2]
        PDict['Number of Complementary Phases'] = y[3]
        PDict['Early Finish1'] = y[4]

        #### FOR NUM OXIDES
        PDict['Number of Oxides (with formation E correction)'] = oxide_no1
        PDict['Oxide List (with formation E correction)'] = oxides_id_withform1
        
        
        x = find_comp(oxide_listdict, compound['unit_cell_formula'], compound['formation_energy_per_atom'], 'Oxide')
        PDict['Complementary Oxide List'] = x[0]
        PDict['Complementary Heat of Oxidation'] = x[1]
        PDict['Lower Formation Energy Than Original Material'] = x[2]
        PDict['Number of Complementary Oxides'] = x[3]
        PDict['Early Finish2'] = x[4]
        
        v_ratio2 = 1000
        for i in x[0]:
            v2 = i['volume']/compound['volume']
            if abs(v2 - 1) < abs(v_ratio2 - 1):
                v_ratio2 = v2
                v_ratio_id2 = i
                
        PDict['Best Volume Ratio'] = v_ratio_id2
        PDict['ID of Best Volume Ratio'] = v_ratio2
    
    return PDict
        


    
if __name__ == '__main__':
        pool = mp.Pool(processes=16)
        print('Calculating Data....')
        DictList = list(tqdm.tqdm(pool.imap(Make_Property_Dict, all_compounds), total=len(all_compounds)))
        
        FinalDF = pd.DataFrame(DictList)

        filename = 'FinalDF_' + str(criteria) + '.pckl'
        f = open(filename, 'wb')
        pickle.dump(FinalDF, f)
        f.close()


        print('Done.')
        
