# The idea of this is to make a list of bam files for every family in my population
# Auther: Darrian Talamantes

# importing modules
import numpy as np

# Pointing to files
Parents_of_Progeny_loc="/home/drt06/Documents/Tall_fescue/Mapping_and_QTL/Mapping_and_QTL/Data/Family_Maker/usable_predicted_parents_double.csv"
Parent_Combos_loc="/home/drt06/Documents/Tall_fescue/Mapping_and_QTL/Mapping_and_QTL/Data/Family_Maker/All_combos.csv"

# reading in files as arrays
Parent_Combos=np.genfromtxt(Parent_Combos_loc, dtype=str, encoding=None, delimiter=",")
Parents_of_Progeny=np.genfromtxt(Parents_of_Progeny_loc, dtype=str, encoding=None, delimiter=",")
Parents_of_Progeny=np.delete(Parents_of_Progeny, 0, 0)

# Directories
Progeny_dir="/scratch/drt83172/Wallace_lab/Mapping_and_QTL/Mapping_and_QTL/Data/Real_Data/Mapped_Reads/"
Parent_dir="/scratch/drt83172/Wallace_lab/Mapping_and_QTL/Mapping_and_QTL/Data/Real_Data/Mapped_Reads_Parents/"

# File extensions
dupped="_dupped.bam"

# loop looks through files and writes out all progeny and parents to a family file.
row_num =     len(Parents_of_Progeny)
for x in Parent_Combos:
    parent1=x[0]
    parent2=x[1]
    parent2=parent2.strip()
    filename=parent1+"x"+parent2+".txt"
    f = open(filename,"w")
    f.write(Parent_dir + parent1 + "\n")
    f.write(Parent_dir + parent2 + "\n")
    for row in range(row_num):
        p1 = Parents_of_Progeny[row][1]
        p2 = Parents_of_Progeny[row][2]
        if (p1 == parent1 or p1 == parent2) and (p2 == parent1 or p2 == parent2):
            f.write(Progeny_dir + Parents_of_Progeny[row][0] + "\n")

    # Create a file called parent1xparent2.txt
    # append all progeny to that file with parent1 and parent 2
    # if col 2 is parent1 or parent2 and col 3 is parent 1 or parent 2 do 
    # write progeny to file 