import os
import csv

def importdic(file):
    d = {}
    with open(file, 'r') as f: 
        csv_reader = csv.reader(f)
        for row in csv_reader:
            key, value = row
            d[key] = value
    return(d)

# Uses a dictionary of file names to rename things. New name on left, old name on right.
def rename_files(dict, data_dir):
    for key, value in dict.items():
        datapath = data_dir + value 
        newname = data_dir + key
        if os.path.exists(datapath):
            os.rename(datapath,newname)


