import pickle
from pymatgen import MPRester
import multiprocessing as mp
import tqdm
import pandas as pd

class StabilityScreener:
    def __init__ (self):
        self.all_compounds = None
        self.properties = None
    
    def load_compounds(self, filename):
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
        self.all_compounds = self.clean_and_check_data(all_compounds_data)
        self.properties = 
    
    def load_compounds(self, api_key, properties, filename=None):
        
        mpr = MPRester(api_key)
        print("Loading Compounds....")
        all_compounds_data = mpr.query({}, properties=properties)
        self.all_compounds = self.clean_and_check_data(all_compounds_data)
        
        if filename is not None:
            file = open(filename, 'wb')
            pickle.dump(all_compounds, file)
            file.close()
    
    def clean_and_check_data(self, data):
        all_compounds = []
        
        print("Cleaning data....")
        for compound in data:
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
    
    def generate_competing_phases(processes, stability_criteria=None, environment=None):
        
        pool = mp.Pool(processes=processes)

        if stability_criteria is None:
            search_space = self.all_compounds
        else:
            print('Finding Stable Phases....')    
            search_space = []
            
            #find all compounds with e_above_hull within specification
            for compound in tqdm.tqdm(self.all_compounds): 
                if abs(compound['e_above_hull']) < stability_criteria/1000:
                    search_space.append(compound)
        
        print('Calculating Data....')
        DictList = list(tqdm.tqdm(pool.imap(Make_Property_Dict, search_space), \
        total=len(search_space)))