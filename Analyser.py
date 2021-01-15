#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  4 18:07:03 2019

@author: NickT
"""

import pickle
import multiprocessing as mp
import tqdm
import pandas as pd


def load_compounds(filename):
    """
    Load Materials Project database saved in .pckl file and check it has
    expected form.
    args:
        filename (str) - name of file containing data
    output:
        list of dicts of data
    """
    print("Loading Compounds....")
    file = open(filename, 'rb')
    all_compounds_data = pickle.load(file)
    all_compounds = []
    print("Cleaning data....")
    for compound in all_compounds_data:
        if None not in compound.values():
            all_compounds.append(compound)

    print("Checking data....")
    discrepancy = 0
    for compound in all_compounds:
        if compound['nsites'] != sum(compound['unit_cell_formula'].values()):
            discrepancy+=1
    assert discrepancy == 0, str(discrepancy) + \
            " compounds have n_sites that don't match their formulas."
    return all_compounds


def find_stable_phases(all_compounds, criteria):
    """
    find all compounds with e_above_hull within given range of zero
    args:
        all_compounds - returned from load_compounds (list of dicts)
        criteria - criteria for a stable phase in meV
    output:
        list of dicts of stable phases
    """
    print('Finding Stable Phases....')    
        
    stable_phase = []
    
    #find all compounds with e_above_hull within 0.05 of 0
    for compound in tqdm.tqdm(all_compounds): 
        if abs(compound['e_above_hull']) < criteria/1000:
            stable_phase.append(compound)
    return stable_phase

def remove_oxide(compounds):
    print('Removing oxides from search....')

    oxide_free = []

    for compound in tqdm.tqdm(compounds):
        if 'O' not in compound['elements']:
            oxide_free.append(compound)
    return oxide_free

######## COMPETING PHASE AND OXIDE CALCULATION ########

def find_comp(stable_oxides, compound_unit_cell, compound_formE, condition, n):
    '''
    Finds complementary oxide or competing phases group and associated total 
    heat of oxidation
    
    args:
        stable_oxides - list of dictionaries of stable oxides or competing 
        phases with lower formation energy than original material
        compound_unit_cell - dict of elements in unit cell of original compound
        ompound_formE - formation energy of original compound
        condition - string dictating whether it is for comp oxide or comp 
        competing phases
        n - forced first list selection (n = 0, 1, 2, 3, ...)
    
    output:
        tuple: (list of dicitionatries of predicted materials, 
        combined formation energy of these materials (with appropriate ratios),
        combined formation energy of these materials (with appropriate ratios) 
        minus original form E,
        number of compounds in complementary group,
        whether the grouping fuction finished early (boolean),
        n - forced first list selection (n = 0, 1, 2, 3, ...))
        
    notes:
        intersect_rank: used to find limiting element by finding ratio of 
        normalised stochiometry between original material and oxide
        
    '''
    result = []
    FinishEarly = False
    #what if positive formE
    
    orig_natoms = sum(compound_unit_cell.values())
    normalised_unit_cell = dict((a, b/orig_natoms) for a, b in \
    compound_unit_cell.items()) #normalise stoichiometry
    
    for oxide in stable_oxides:
        oxide['el_weight'] = dict((a, b/oxide['nsites']) for a, b in \
        oxide['unit_cell_formula'].items()) #normalise stoichiometry
        
        if condition == 'Oxide':
            del oxide['el_weight']['O']
        
        oxide['ranker'] = dict((a, b/normalised_unit_cell[a]) for a, b in \
        oxide['el_weight'].items()) #find greedy ranking parameter
        
        # define how much material is used in total per unit cell
        oxide['ranking_no'] = sum(oxide['ranker'].values())
    
    # Order by energy per unit used up
    sort_oxides = sorted(stable_oxides, key = lambda oxide: \
    (oxide['formation_energy_per_atom']/oxide['ranking_no']))
        
    sort_oxides1 = sort_oxides[:]
    
    total_formE = 0
    counter = 0
    #if all atoms in unit cell not yet accounted for
    while sum(normalised_unit_cell.values()) != 0 and sort_oxides1 != []:         
        if counter == 0:    
            oxide = sort_oxides1[n]
        else: 
            oxide = sort_oxides1[0] #to allow forced initial choice
        
        # elements in oxide and orig. material
        intersection = list(set(oxide['elements']).\
        intersection(normalised_unit_cell.keys()))
        # shouldnt we remove O from intersection???
        if len(intersection) == 0:
            print(compound_unit_cell)
            print(oxide['unit_cell_formula'])
            print(oxide['nsites'])
        intersect_rank = {}

        for element in intersection:
            # same as 1/ranker values
            intersect_rank[element] = normalised_unit_cell[element]/ \
            (oxide['unit_cell_formula'][element]/oxide['nsites'])
        
        #find limiting element   
        limiting_element = min(intersect_rank, key=intersect_rank.get) 
        ratio = intersect_rank[limiting_element] #(value)
        oxide['ratio'] = ratio # For PBR calculation
        used_up_elements = []
        for element in intersection:
            
            normalised_unit_cell[element] = normalised_unit_cell[element] - \
            (ratio * oxide['unit_cell_formula'][element]/oxide['nsites'])
            
            #inequality because of != 0 problem
            if abs(normalised_unit_cell[element]) < 0.0001: 
                used_up_elements.append(element)
                
        result.append(oxide)
        sort_oxides1.remove(oxide)
        total_formE += oxide['formation_energy_per_atom']*ratio
        
        #remove oxides in list which arent useful (dont have new elements)
        sort_oxides1 = [oxide for oxide in sort_oxides1 if \
        len(set(oxide['elements']).intersection(used_up_elements)) == 0]
        counter += 1
        
    #inequality because of != 0 problem
    if len(sort_oxides1) == 0 and \
                        abs(sum(normalised_unit_cell.values())) > 0.0001: 
        
        FinishEarly = True
                
    return (result, total_formE, total_formE-compound_formE, \
                len(result), FinishEarly, n)


def forced_choice(stable_oxides, compound_unit_cell, compound_formE, condition):
    '''
    Function which forces the find_comp function to adjust its first choice 
    in the oxide list.
    
    When there are 3 or more oxides find_comp will be run 3 times, forcing 
    1st, 2nd and 3rd on the ordered list to be selected first. The function 
    then result with the minimum total_formE-compound_formE.
    
    When there are 1 or 2 oxides find_comp will only be called once or twice, 
    respectively.
    
    args:
        stable_oxides - list of dictionaries of stable oxides or competing 
        phases with lower formation energy than original material
        compound_unit_cell - dict of elements in unit cell of original compound
        ompound_formE - formation energy of original compound
        condition - string dictating whether it is for comp oxide or comp 
        competing phases
        
    output:
        tuple: (list of dicitionatries of predicted materials, 
        combined formation energy of these materials (with appropriate ratios),
        combined formation energy of these materials (with appropriate ratios) 
        minus original form E,
        number of compounds in complementary group,
        whether the grouping fuction finished early (boolean),
        index of oxide first selected by find_comp (n) which gave this result)
        
    '''
    length = len(stable_oxides)
    
    if length == 1:
        ans = find_comp(stable_oxides, compound_unit_cell, compound_formE, \
        condition, 0)
        
    elif length == 2:
        x1 = find_comp(stable_oxides, compound_unit_cell, compound_formE, \
        condition, 0)
        x2 = find_comp(stable_oxides, compound_unit_cell, compound_formE, \
        condition, 1)
        
        ans = min([x1, x2], key = lambda x: x[2])
    
    else:
        x1 = find_comp(stable_oxides, compound_unit_cell, compound_formE, \
        condition, 0)
        x2 = find_comp(stable_oxides, compound_unit_cell, compound_formE, \
        condition, 1)
        x3 = find_comp(stable_oxides, compound_unit_cell, compound_formE, \
        condition, 2)
        
        ans = min([x1, x2, x3], key = lambda x: x[2])
    
    return ans

        



    #### FOR TESTING FIND_OXIDES
def find_oxides_test():
            
    ABCO4 = {'elements': ['A', 'B', 'C', 'O'], 'formation_energy_per_atom':\
    -750, 'nsites':7, 'unit_cell_formula':{'A':1, 'B':1, 'C':1, 'O':4}}
    
    AO = {'elements': ['A', 'O'], 'formation_energy_per_atom': -100, \
    'nsites':8, 'unit_cell_formula':{'A':4, 'O':4}}
    
    BO2 = {'elements': ['B', 'O'], 'formation_energy_per_atom': -100, \
    'nsites':6, 'unit_cell_formula':{'B':2, 'O':4}}
    
    C2O = {'elements': ['C', 'O'], 'formation_energy_per_atom': -300, \
    'nsites':24, 'unit_cell_formula':{'C':16, 'O':8}}
    
    A2BO6 = {'elements': ['A', 'B', 'O'], 'formation_energy_per_atom': -380, \
    'nsites':9, 'unit_cell_formula':{'A':2, 'B':1, 'O':6}}
    
    A2CO4 = {'elements': ['A', 'C', 'O'], 'formation_energy_per_atom': -620, \
    'nsites':63, 'unit_cell_formula':{'A':18, 'C':9, 'O':36}}
    
    original = {'A':4, 'B':8, 'C':100}
   
    oxides = [ABCO4, AO, BO2, C2O, A2BO6, A2CO4]
    # Expected result for first pick is dH = -78.9286 eV/atom
    print(find_comp(oxides, original, -400, 'Oxide', 0))
    print("***********")
    print(forced_choice(oxides, original, -400, 'Oxide'))    
                        
      

def Make_Property_Dict(compound):
    '''
    Function to be iterated over all compounds.
    '''
    PDict = {}
    
    if abs(compound['e_above_hull']) < criteria/1000: #if stable 
        
        #### FOR NUM PHASES
        competing_phases_id_withform = []
        competing_phase_no = 0
        comp_listdict =[]
        
        
        #### FOR NUM OXIDES
        oxide_no = 0
        oxides_id_withform = []
        oxide_listdict = []
            
        elements = compound['elements']
        
        for i in stable_phase:
            #### FOR NUM PHASES
            if set(i['elements']).issubset(elements):
                comp_listdict.append(i) #for find_comp
                
                if i['formation_energy_per_atom'] < \
                                compound['formation_energy_per_atom']:
                    #find all other phases containing just those elements                
                    competing_phase_no +=1
                    competing_phases_id_withform.append(i['task_id'])
                    
                    
            #### FOR NUM OXIDES
            if 'O' in i['elements']:
                el = i['elements'][:]
                el.remove('O')
                #o_sites = i['unit_cell_formula']['O']
                if set(el).issubset(elements) and len(el) != 0:
                    oxide_listdict.append(i) #for find_comp
                    
                    if i['formation_energy_per_atom'] < \
                                    compound['formation_energy_per_atom']:
                        oxide_no += 1
                        oxides_id_withform.append(i['task_id'])
        
        
        
        #### FOR NUM PHASES
        
        PDict['task_id'] = compound['task_id']
        PDict['Formula'] = compound['pretty_formula']
        assert type(PDict['Formula']) == str
        PDict['Bandgap /eV'] = compound['band_gap']
        PDict['e_above_hull'] = compound['e_above_hull']
        PDict['elements'] = compound['elements']
        PDict['Competing Phase List (with formation E correction)'] = \
        competing_phases_id_withform
        
        y = forced_choice(comp_listdict, compound['unit_cell_formula'], \
        compound['formation_energy_per_atom'], 'NotOx')
        
        PDict['Complementary Competing Phase List'] = y[0]
        PDict['Heat of Decomposition'] = y[2]
        PDict['Early Finish1'] = y[4]
        PDict['Index of First Phase Selected'] = y[5]
        
        #### FOR NUM OXIDES
        PDict['Number of Oxides (with formation E correction)'] = oxide_no
        PDict['Oxide List (with formation E correction)'] = oxides_id_withform
        
        x = forced_choice(oxide_listdict, compound['unit_cell_formula'], \
        compound['formation_energy_per_atom'], 'Oxide')
        
        PDict['Complementary Oxide List'] = x[0]
        PDict['Heat of Oxidation'] = x[2]
        PDict['Early Finish2'] = x[4]
        PDict['Index of First Oxide Selected'] = x[5]


        vol_ratio = []
        
        for i in x[0]:
            vol_ratio.append((i['volume']*i['ratio']*compound['nsites']) / \
            (compound['volume']*i['nsites'])) # Assuming good diffusion
        # Use to see whether top two oxides have PBR 1-2        
        PDict['Volume Ratios'] = vol_ratio

    
    return PDict
        


    
if __name__ == '__main__':
    
    all_compounds = load_compounds("MPDatabase.pckl")
    criteria = 50 # criteria for stable phases in meV
    stable_phase = find_stable_phases(all_compounds, criteria)
    
    print(len(stable_phase))
    
    oxide_free = remove_oxide(stable_phase)
    print(oxide_free)

    pool = mp.Pool(processes=32)
    print('Calculating Data....')
    DictList = list(tqdm.tqdm(pool.imap(Make_Property_Dict, oxide_free), \
    total=len(stable_phase)))
    
    FinalDF = pd.DataFrame(DictList)
    
    print('Saving result....')
    filename = 'FinalDF_noox_' + str(criteria) + '.pckl'
    f = open(filename, 'wb')
    pickle.dump(FinalDF, f)
    f.close()
    '''
    #find_oxides_test()
    all_compounds = load_compounds("MPDatabase.pckl")
    for i in all_compounds:
        if i["pretty_formula"] == "EuAg":
            print(i)
    #for i in stable_phase:
        #### FOR NUM PHASES
    #    if set(i['elements']).issubset(elements):
    #        comp_listdict.append(i) #for find_comp
    '''     
    print('Done.')
        
