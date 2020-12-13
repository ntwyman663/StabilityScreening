import pickle
import pandas as pd

print("Loading Compounds....")
file = open('MPDatabase.pckl', 'rb')
database = pickle.load(file)

print(type(database))
print("Number of Rows: ", len(database.index))

