import pickle
import pandas as pd

print("Loading Compounds....")
file = open('MPDatabase.pckl', 'rb')
database = pickle.load(file)

print(type(database))
print("Number of Rows: ", len(database))

attribute_counts = [0 for i in range(100)]
for i in database:
    attribute_counts[len(i.values())] += 1

count_dict = {}
for i, att in enumerate(attribute_counts):
    if att!=0:
        print(att, " compounds with ", i, " attributes")

discrepancy = 0
for compound in database:
    if compound['nsites'] != sum(compound['unit_cell_formula'].values()):
        discrepancy+=1
assert discrepancy == 0, str(discrepancy) + \
        " compounds have n_sites that don't match their formulas."
count = 0
for i in database:
    if type(i['pretty_formula']) != str:
        count+=1
print(count)


file = open('FinalDictList_50.pckl', 'rb')
DL = pickle.load(file)
DF = pd.DataFrame(DL)

print("Saving...")
file = open('FinalDF_50.pckl', 'wb')
pickle.dump(DF, file)
file.close()
print("Done.")
