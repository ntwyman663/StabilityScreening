"""
Created on Mon Jan 28 23:08:26 2019

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
                          'elements', 'volume', 'formation_energy_per_atom', 'band_gap'])

criteria = float(input("Enter Stable Phase Criteria in meV: ")) 

def find_stable_phases(compound):
    '''
    find all compounds with e_above_hull within given range of zero
    '''
    if abs(compound['e_above_hull']) < criteria/1000:
        return compound


print('Finding Stable Phases....')    
    
stable_phase = []

for compound in tqdm.tqdm(all_compounds): #find all compounds with e_above_hull within 0.05 of 0
    if abs(compound['e_above_hull']) < criteria/1000:
        stable_phase.append(compound)

#pool = mp.Pool(processes=1)
#
#stable_phase = list(tqdm.tqdm(pool.imap(find_stable_phases, all_compounds), total=86680))


######## COMPETING PHASE AND OXIDE CALCULATION ########

def find_oxides(stable_oxides, compound_elements, compound_formE):
    '''
    Finds the minimum number stable oxides which need to form to account
    for every element in compound element.
    
    e.g. for compound ABC if list of stable oxides ordered with lowest 
    formation energy per element first is [AO, ABO, BO, BCO, CO]
    the result from this function would be [AO, ABO, BCO]. This covers all elements
    
    Oxides weighted by formE/ no. of elements 
    
    args:
        stable_oxides = list of dictionaries of stable oxides with lower formation
        energy than original material
        compound_elements = list of elements (strings) which are in original compound
        **compound_formE = formation energy of compound**
    
    output:
        tuple: (list of dicitionatries of predicted oxides, 
        combined formation energy of produced oxides,
        whether this combined formE is lower than that of original material (boolean))
        
    '''
    result = []
    #what if positive formE
    sort_oxides = sorted(stable_oxides, key = lambda oxide: oxide['formation_energy_per_atom']/
                         (len(oxide['elements'])-1))
    
    compound_elements1 = compound_elements[:]
    total_formE = 0
    for oxide in sort_oxides:
        if compound_elements1 != []:
            ox_elements = oxide['elements']
            x = 0
            for i in ox_elements:
                if i in compound_elements1:
                    compound_elements1.remove(i)
                    x += 1
            if x >= 1:
                result.append(oxide)
                total_formE += oxide['formation_energy_per_atom']
                
        else: break
                
    return (result, total_formE, total_formE<compound_formE)

    #### FOR TESTING FIND_OXIDES
        
#ABCO = {'elements': ['A', 'B', 'C', 'O'], 'formation_energy_per_atom': -200}
#AO = {'elements': ['A', 'O'], 'formation_energy_per_atom': -100}
#BO = {'elements': ['B', 'O'], 'formation_energy_per_atom': -100}
#CO = {'elements': ['C', 'O'], 'formation_energy_per_atom': -100}
#
#listt = [ABCO, AO, BO, CO]
                
#find_oxides(listt, ['A', 'B', 'C'], -200)    
                    
    ####  
      

def Make_Property_Dict(compound):
    '''
    Function to be iterated over all compounds.
    '''
    PDict = {}
    
    global stable_phase
    
    if abs(compound['e_above_hull']) < criteria/1000: #if stable 
        
        #### FOR NUM PHASES
        competing_phase_id1 = []
        competing_phases_id_withform1 = []
        competing_phase_no = 0
        competing_phase_no1 = 0
        
        
        #### FOR NUM OXIDES
        v_ratio = 0
        oxides_id1 = []
        oxide_no = 0
        v_ratio_id1 = 'n/a'
        
        v_ratio2 = 0
        oxide_no1 = 0
        oxides_id_withform1 = []
        v_ratio_id2 = 'n/a'
        oxide_listdict = []
            
        elements = compound['elements']
        
        for i in stable_phase:
            #### FOR NUM PHASES
            if set(i['elements']).issubset(elements): #find all other phases containing just those elements
                competing_phase_no += 1 #add number
                competing_phase_id1.append(i['task_id']) #add phase ID
                
                #make seperate list including formation energy criteria
                if i['formation_energy_per_atom'] < compound['formation_energy_per_atom']:
                    competing_phase_no1 +=1
                    competing_phases_id_withform1.append(i['task_id'])
                    
                    
            #### FOR NUM OXIDES
            if 'O' in i['elements']:
                el = i['elements'][:]
                el.remove('O')
                if set(el).issubset(elements) and el != []:
                    v1 = i['volume']/compound['volume']
                    oxide_no += 1
                    oxides_id1.append(i['task_id'])
                    oxide_listdict.append(i)
    
                    if abs(v1 - 1) < abs(v_ratio - 1):
                        v_ratio = v1
                        v_ratio_id1 = i['task_id']
                    
                    if i['formation_energy_per_atom'] < compound['formation_energy_per_atom']:
                        oxide_no1 += 1
                        v2 = i['volume']/compound['volume']
                        oxides_id_withform1.append(i['task_id'])
                        
                        if abs(v2 - 1) < abs(v_ratio2 - 1):
                            v_ratio2 = v2
                            v_ratio_id2 = i['task_id']
        
        
        
        #### FOR NUM PHASES
        #### remove entry of itself from competing phases ###
        competing_phase_no -= 1
        competing_phase_id1.remove(compound['task_id'])        
        ##
        
        PDict['task_id'] = compound['task_id']
        PDict['Formula'] = compound['pretty_formula']
        PDict['Bandgap /eV'] = compound['band_gap']
        PDict['Number of Competing Phases'] = competing_phase_no
        PDict['Competing Phase List'] = competing_phase_id1
        
        PDict['Competing Phase Number (with formation E correction)'] = competing_phase_no1
        PDict['Competing Phase List (with formation E correction)'] = competing_phases_id_withform1
        
        
        #### FOR NUM OXIDES
        PDict['Number of Oxides'] = oxide_no
        PDict['Oxide List'] = oxides_id1
        PDict['Material:Oxide Volume Ratio Closest to 1'] = v_ratio
        PDict['Oxide with Best Ratio'] = v_ratio_id1
        
        PDict['Number of Oxides (with formation E correction)'] = oxide_no1
        PDict['Oxide List (with formation E correction)'] = oxides_id_withform1
        PDict['Material:Oxide Ratio Closest to 1 (with formation E correction)'] = v_ratio_id2
        PDict['Oxide with Best Ratio (with formation E correction)'] = v_ratio2
        
        x = find_oxides(oxide_listdict, elements, compound['formation_energy_per_atom'])
        PDict['Complementary Oxide List'] = x[0]
        PDict['Lower Formation Energy Than Original Material'] = x[2]

        
    else:
        #### FOR NUM PHASES
        PDict['task_id'] = compound['task_id']
        PDict['Formula'] = compound['pretty_formula']
        PDict['Bandgap /eV'] = compound['band_gap']
        PDict['Number of Competing Phases'] = 'unstable'
        PDict['Competing Phase List'] = 'unstable'
        
        PDict['Competing Phase Number (with formation E correction)'] = 'unstable'
        PDict['Competing Phase List (with formation E correction)'] = 'unstable'
        
        
        #### FOR NUM OXIDES
        PDict['Number of Oxides'] = 'unstable'
        PDict['Oxide List'] = 'unstable'
        PDict['Material:Oxide Volume Ratio Closest to 1'] = 'unstable'
        PDict['Oxide with Best Ratio'] = 'unstable'
        
        PDict['Number of Oxides (with formation E correction)'] = 'unstable'
        PDict['Oxide List (with formation E correction)'] = 'unstable'
        PDict['Material:Oxide Ratio Closest to 1 (with formation E correction)'] = 'unstable'
        PDict['Oxide with Best Ratio (with formation E correction)'] = 'unstable'
        

        PDict['Complementary Oxide List'] = 'unstable'
        PDict['Lower Formation Energy Than Original Material'] = 'unstable'
    
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

        #f = open('stable_phase_50meV.pckl', 'wb')
        #pickle.dump(stable_phase, f)
        #f.close()	

        print('Done.')
