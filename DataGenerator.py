import pickle
from pymatgen import MPRester

mat_api_key = 'JWV6Fi4f6VfxROtHO2uP'

mpr = MPRester(mat_api_key)

print("Loading Compounds....")
all_compounds = mpr.query({}, properties=["task_id", "pretty_formula", 'e_above_hull',
                          'elements', 'volume', 'formation_energy_per_atom', 'band_gap',
                           'nsites', 'unit_cell_formula'])

file = open('MPDatabase.pckl', 'wb')
pickle.dump(all_compounds, file)
file.close()
