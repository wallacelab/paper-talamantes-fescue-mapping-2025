#! Python
# Objective is to load in a vcf file as a data frame and create tables that can load into R for visulization

# test file
# /home/drt06/Documents/Tall_fescue/Mapping_and_QTL/Mapping_and_QTL/Data/Real_Data/VCF/haplotypes/beagle_315x320.vcf
# /home/drt06/Documents/Tall_fescue/Mapping_and_QTL/Mapping_and_QTL/Data/Real_Data/Lists/usable_predicted_parents_double.csv
# /home/drt06/Documents/Tall_fescue/Mapping_and_QTL/Mapping_and_QTL/Data/Real_Data/VCF/haplotypes/Haplotype_Eval_Data

import pandas as pd
import numpy as np
import argparse
import io
import os


def main():
    # vaeiables 
    data_folder = "/home/drt06/Documents/Tall_fescue/Mapping_and_QTL/Mapping_and_QTL/Data/Real_Data/VCF/haplotypes/Haplotype_Eval_Data"
    og_vcf_loc = "/home/drt06/Documents/Tall_fescue/Mapping_and_QTL/Mapping_and_QTL/Data/Real_Data/VCF/all_maf_05_min_65_no_indels_not_multiallelic_315x320.vcf"
    full_vcf_loc = "/home/drt06/Documents/Tall_fescue/Mapping_and_QTL/Mapping_and_QTL/Data/Real_Data/VCF/all_maf_05_min_65_no_indels_not_multiallelic.vcf"
    # Creating variables
    window_size = 15
    step = 3

    # Loading in files
    # args.vcf is the vcf file
    # args.parent is a file of predicted parents
    args = parse_args()
    vcf_file = VCF_importer(args.vcf)
    parent_file = csv_importer(args.parent) 
    og_vcf = VCF_importer(og_vcf_loc)
    full_vcf = VCF_importer(full_vcf_loc)
    # haplotype_percents(vcf_file, parent_file, data_folder, window_size, step)
    # haplotype_percents_parents(vcf_file, parent_file, data_folder, window_size, step)
    # allelfrequncie(vcf_file,data_folder)
    # vcf_vs_beagle(og_vcf, vcf_file, data_folder)
    find_geno_depth(full_vcf, data_folder)

# This takes the depth feild and divides it by number of progeny without missing data.
def find_geno_depth(full_vcf, data_folder):
    results_array = np.zeros((len(full_vcf), 2), dtype='<U30')
    full_vcf['POS'] = full_vcf['POS'].astype(str)
    full_vcf['ID'] = pd.concat([full_vcf['CHROM'], full_vcf['POS']], axis=1).apply(lambda row: ''.join(row), axis=1)
    for row in range(0, len(full_vcf)):
        ID = full_vcf.iloc[row,2]
        results_array[row,0] = ID
        # Extracting the depth from VCF file
        infoline = full_vcf.iloc[row, 7]
        infoline = infoline.split(";")
        dpfeild = infoline[0]
        dpfeild = int(dpfeild[3:])
        not_nullcount = 0
        
        for col in range(9,len(full_vcf.columns)):
            genotype = full_vcf.iloc[row,col]
            genotype = genotype[:3]
            if (genotype!="./."):
                # print(dpfeild, genotype)
                not_nullcount = not_nullcount + 1
        results_array[row,1] = dpfeild/not_nullcount
    file_name = "genotype_depth.csv"    
    file_name = os.path.join(data_folder, file_name)
    numeric_indices = np.char.isnumeric(results_array)
    # Convert the numeric values to float
    results_array[numeric_indices] = results_array[numeric_indices].astype(float)
    # Save the array with appropriate format specifier
    np.savetxt('output.txt', results_array, fmt='%s')
    np.savetxt(file_name, results_array, delimiter=",", fmt='%s')

def pseudo_testcross_detector(full_vcf, data_folder):
    print("Hi")


# this function will compare the output of beagle to its regular vcf counterpart to see how much it changes.
def vcf_vs_beagle(og_vcf, beagle_vcf, data_folder):
    # Creating unique identifiers
    og_vcf['POS'] = og_vcf['POS'].astype(str)
    og_vcf['ID'] = pd.concat([og_vcf['CHROM'], og_vcf['POS']], axis=1).apply(lambda row: ''.join(row), axis=1)
    beagle_vcf['POS'] = beagle_vcf['POS'].astype(str)
    beagle_vcf['ID'] = pd.concat([beagle_vcf['CHROM'], beagle_vcf['POS']], axis=1).apply(lambda row: ''.join(row), axis=1)
    identifier = beagle_vcf[["ID"]]
    og_vcf2 = pd.merge(og_vcf, identifier, on='ID', how='inner')
    
    print(og_vcf2.head())
    print(og_vcf2.shape)
    print(beagle_vcf.shape)
    results_array = np.zeros((len(beagle_vcf), 6))
    results_array = results_array.astype(object)

    # sets up a for loop to go through both files
    for row  in range(0, len(beagle_vcf)):
        same = 0
        Homo2Het = 0 
        Het2Homo = 0
        Het2Het = 0
        filled = 0
        for col in range(9,len(beagle_vcf.columns)):
            bealge = beagle_vcf.iloc[row,col]
            og = og_vcf2.iloc[row,col]
            og = og[:3]
            og = og.replace("/","|")
            print(og," ", bealge)

            # Checking if 
            ogID = og_vcf2.iloc[row,2]
            beagleID = beagle_vcf.iloc[row,2]

            # This massive loop checks if the IDs are the same then places it in the results table
            # It then checks for 3 different outcomes and adds one to them wherever it hits
            if (ogID == beagleID):
                if (bealge == og):
                    print("Its the same", same)
                    same = same + 1
                    results_array[row,1] = same

                elif og == ".|."  and bealge != ".|.":
                    print("Thats a 0", filled)
                    filled = filled + 1
                    results_array[row,2] = filled

                elif (og == "1|1" or og == "0|0") and (bealge == "1|0" or bealge == "0|1"):
                    print("Thats homo!", Homo2Het)
                    Homo2Het = Homo2Het + 1
                    results_array[row,3] = Homo2Het

                elif (og == "0|1" or og == "1|0") and (bealge == "1|1" or bealge == "0|0"):
                    print("Thats homo!", Het2Homo)
                    Het2Homo = Het2Homo + 1
                    results_array[row,4] = Het2Homo 

                elif (og == "0|1" and bealge == "1|0" ):
                    print("Thats hetero!", Het2Het)
                    Het2Het = Het2Het + 1
                    results_array[row,5] = Het2Het

                elif (og == "1|0" and bealge == "0|1" ):
                    print("Thats hetero!", Het2Het)
                    Het2Het = Het2Het + 1
                    results_array[row,5] = Het2Het
                else:
                    y =+ 1

    resultsdf = pd.DataFrame(results_array, columns =['ID','Same', 'Filled', 'Homo2Het', 'Het2Homo', 'Het2Het'])
    file_name = "beagle_vs_vcf.csv"
    file_name = os.path.join(data_folder, file_name)
    resultsdf.to_csv(file_name)


 


# This creates a file with percent similarity for each sliding window to each progeny to parental haplotype and saves the data.
def haplotype_percents(vcf_file, parent_file, data_folder, window_size, step):

    # Define Variables
    plant_names = vcf_file.columns.tolist()
    tot_prog = len(parent_file.index)

    # This loop iterates through the parent file to create variables for the current working individuals
    for prog_num in range (0, 3):
        current_proj = parent_file.iloc[prog_num]["Progeny"]
        current_parent1 = parent_file.iloc[prog_num]["Parent1"]
        current_parent2 = parent_file.iloc[prog_num]["Parent2"]
    
    test_prog = "315-4-43"
    test_parent1 = "320"
    test_parent2 = "315"
    # This loop checks if progeny profile is present in the vcf then starts to work with it and its parents.
    if test_prog in plant_names:
        print(f"{test_prog} is in the list.")
        result_df = pd.DataFrame(columns=[ 'Chrom','ProgAxParent1A', 'ProgAxParent1B', 'ProgAxParent2A', 'ProgAxParent2B', 'ProgBxParent1A', 'ProgBxParent1B', 'ProgBxParent2A', 'ProgBxParent2B'])
        percent_file = os.path.join(data_folder, test_prog)
        x=0
        Chrome_list = vcf_file['CHROM'].sort_values().unique()
        for chrom in range(0, len(Chrome_list)):
            current_chrome=Chrome_list[chrom]
            vcf_file_subset = vcf_file[vcf_file['CHROM'] == current_chrome]
            if len(vcf_file_subset) >= (window_size + step*2):
                
                for i in range(0, len(vcf_file_subset) - window_size + 1, step):
                    # Extract windows from each column
                    window_prog = vcf_file_subset[test_prog].iloc[i:i+window_size]
                    window_parent1 = vcf_file_subset[test_parent1].iloc[i:i+window_size]
                    window_parent2 = vcf_file_subset[test_parent2].iloc[i:i+window_size]

                    # splitting windows
                    prog_split_data = window_prog.str.split("|", expand=True)
                    prog_split_data.columns = ['ProgA', 'ProgB']
                    parent1_split_data = window_parent1.str.split("|", expand=True)
                    parent1_split_data.columns = ['Parent1A', 'Parent1B']
                    parent2_split_data = window_parent2.str.split("|", expand=True)
                    parent2_split_data.columns = ['Parent2A', 'Parent2B']

                    """if x == 7:
                        print("progeny ", prog_split_data)
                        print("parent1 ", parent1_split_data)
                        print("parent1 ", parent2_split_data)
                    x= x+1 """

                    # Doing comparisons between sliding windows and saving the results to a data frame 
                    # Count the number of matching positions
                    matching_positions_ProgAxParent1A = (prog_split_data['ProgA'] == parent1_split_data['Parent1A']).sum()
                    matching_positions_ProgAxParent1B = (prog_split_data['ProgA'] == parent1_split_data['Parent1B']).sum()
                    matching_positions_ProgAxParent2A = (prog_split_data['ProgA'] == parent2_split_data['Parent2A']).sum()
                    matching_positions_ProgAxParent2B = (prog_split_data['ProgA'] == parent2_split_data['Parent2B']).sum()

                    matching_positions_ProgBxParent1A = (prog_split_data['ProgB'] == parent1_split_data['Parent1A']).sum()
                    matching_positions_ProgBxParent1B = (prog_split_data['ProgB'] == parent1_split_data['Parent1B']).sum()
                    matching_positions_ProgBxParent2A = (prog_split_data['ProgB'] == parent2_split_data['Parent2A']).sum()
                    matching_positions_ProgBxParent2B = (prog_split_data['ProgB'] == parent2_split_data['Parent2B']).sum()

                    # Count the total number of positions
                    total_positions = len(prog_split_data)
                    # Calculate the percentage similarity
                    percentage_ProgAxParent1A = (matching_positions_ProgAxParent1A / total_positions) 
                    percentage_ProgAxParent1B = (matching_positions_ProgAxParent1B / total_positions) 
                    percentage_ProgAxParent2A = (matching_positions_ProgAxParent2A / total_positions) 
                    percentage_ProgAxParent2B = (matching_positions_ProgAxParent2B / total_positions) 

                    percentage_ProgBxParent1A = (matching_positions_ProgBxParent1A / total_positions) 
                    percentage_ProgBxParent1B = (matching_positions_ProgBxParent1B / total_positions) 
                    percentage_ProgBxParent2A = (matching_positions_ProgBxParent2A / total_positions) 
                    percentage_ProgBxParent2B = (matching_positions_ProgBxParent2B / total_positions) 
                    
                    new_data = pd.DataFrame({
                        'Chrom': [current_chrome],
                        'ProgAxParent1A': [round(percentage_ProgAxParent1A,2)],
                        'ProgAxParent1B': [round(percentage_ProgAxParent1B,2)],
                        'ProgAxParent2A': [round(percentage_ProgAxParent2A,2)],
                        'ProgAxParent2B': [round(percentage_ProgAxParent2B,2)],

                        'ProgBxParent1A': [round(percentage_ProgBxParent1A,2)],
                        'ProgBxParent1B': [round(percentage_ProgBxParent1B,2)],
                        'ProgBxParent2A': [round(percentage_ProgBxParent2A,2)],
                        'ProgBxParent2B': [round(percentage_ProgBxParent2B,2)]
                    })
                    # Concatenate the existing DataFrame with the new data
                    result_df = pd.concat([result_df, new_data], ignore_index=True)

                    result_df.to_csv(percent_file)
            else:
                print(current_chrome, " does not have enough snps")       
    else:
        print("progeny is missing") 



def haplotype_percents_parents(vcf_file, parent_file, data_folder, window_size, step):
    # Define Variables
    plant_names = vcf_file.columns.tolist()
    tot_prog = len(parent_file.index)

    # This loop iterates through the parent file to create variables for the current working individuals
    for prog_num in range (0, 3):
        current_proj = parent_file.iloc[prog_num]["Progeny"]
        current_parent1 = parent_file.iloc[prog_num]["Parent1"]
        current_parent2 = parent_file.iloc[prog_num]["Parent2"]
    
    test_prog = "315-4-43"
    test_parent1 = "320"
    test_parent2 = "315"
    cross = "315x320"
    # This loop checks if progeny profile is present in the vcf then starts to work with it and its parents.
    if test_prog in plant_names:
        print(f"{test_prog} is in the list.")
        result_df = pd.DataFrame(columns=[ 'Chrom','Parent1AxParent1A', 'Parent1AxParent1B', 'Parent1AxParent2A', 'Parent1AxParent2B', 'Parent1BxParent1A', 'Parent1BxParent1B', 'Parent1BxParent2A', 'Parent1BxParent2B'])
        percent_file = os.path.join(data_folder, cross)
        x=0
        Chrome_list = vcf_file['CHROM'].sort_values().unique()
        for chrom in range(0, len(Chrome_list)):
            current_chrome=Chrome_list[chrom]
            vcf_file_subset = vcf_file[vcf_file['CHROM'] == current_chrome]
            if len(vcf_file_subset) >= (window_size + step*2):
                
                for i in range(0, len(vcf_file_subset) - window_size + 1, step):
                    # Extract windows from each column
                    window_prog = vcf_file_subset[test_prog].iloc[i:i+window_size]
                    window_parent1 = vcf_file_subset[test_parent1].iloc[i:i+window_size]
                    window_parent2 = vcf_file_subset[test_parent2].iloc[i:i+window_size]

                    # splitting windows
                    prog_split_data = window_prog.str.split("|", expand=True)
                    prog_split_data.columns = ['ProgA', 'ProgB']
                    parent1_split_data = window_parent1.str.split("|", expand=True)
                    parent1_split_data.columns = ['Parent1A', 'Parent1B']
                    parent2_split_data = window_parent2.str.split("|", expand=True)
                    parent2_split_data.columns = ['Parent2A', 'Parent2B']

                    # Doing comparisons between sliding windows and saving the results to a data frame 
                    # Count the number of matching positions
                    matching_positions_ProgAxParent1A = (parent1_split_data['Parent1A'] == parent1_split_data['Parent1A']).sum()
                    matching_positions_ProgAxParent1B = (parent1_split_data['Parent1A'] == parent1_split_data['Parent1B']).sum()
                    matching_positions_ProgAxParent2A = (parent1_split_data['Parent1A'] == parent2_split_data['Parent2A']).sum()
                    matching_positions_ProgAxParent2B = (parent1_split_data['Parent1A'] == parent2_split_data['Parent2B']).sum()

                    matching_positions_ProgBxParent1A = (parent1_split_data['Parent1B'] == parent1_split_data['Parent1A']).sum()
                    matching_positions_ProgBxParent1B = (parent1_split_data['Parent1B'] == parent1_split_data['Parent1B']).sum()
                    matching_positions_ProgBxParent2A = (parent1_split_data['Parent1B'] == parent2_split_data['Parent2A']).sum()
                    matching_positions_ProgBxParent2B = (parent1_split_data['Parent1B'] == parent2_split_data['Parent2B']).sum()

                    # Count the total number of positions
                    total_positions = len(prog_split_data)
                    # Calculate the percentage similarity
                    percentage_Parent1AxParent1A = (matching_positions_ProgAxParent1A / total_positions) 
                    percentage_Parent1AxParent1B = (matching_positions_ProgAxParent1B / total_positions) 
                    percentage_Parent1AxParent2A = (matching_positions_ProgAxParent2A / total_positions) 
                    percentage_Parent1AxParent2B = (matching_positions_ProgAxParent2B / total_positions) 

                    percentage_Parent1BxParent1A = (matching_positions_ProgBxParent1A / total_positions) 
                    percentage_Parent1BxParent1B = (matching_positions_ProgBxParent1B / total_positions) 
                    percentage_Parent1BxParent2A = (matching_positions_ProgBxParent2A / total_positions) 
                    percentage_Parent1BxParent2B = (matching_positions_ProgBxParent2B / total_positions) 
                    
                    new_data = pd.DataFrame({
                        'Chrom': [current_chrome],
                        'Parent1AxParent1A': [round(percentage_Parent1AxParent1A,2)],
                        'Parent1AxParent1B': [round(percentage_Parent1AxParent1B,2)],
                        'Parent1AxParent2A': [round(percentage_Parent1AxParent2A,2)],
                        'Parent1AxParent2B': [round(percentage_Parent1AxParent2B,2)],

                        'Parent1BxParent1A': [round(percentage_Parent1BxParent1A,2)],
                        'Parent1BxParent1B': [round(percentage_Parent1BxParent1B,2)],
                        'Parent1BxParent2A': [round(percentage_Parent1BxParent2A,2)],
                        'Parent1BxParent2B': [round(percentage_Parent1BxParent2B,2)]
                    })
                    # Concatenate the existing DataFrame with the new data
                    result_df = pd.concat([result_df, new_data], ignore_index=True)

                    result_df.to_csv(percent_file)
            else:
                print(current_chrome, " does not have enough snps")       
    else:
        print("progeny is missing")


def allelfrequncie(vcf_file, data_folder):
    filename = "allel_frequencies.csv"
    freq_file = os.path.join(data_folder,filename)
    result_df = pd.DataFrame(columns=[ 'Site'])
    vcf_file[["CHROM", "POS"]] = vcf_file[["CHROM", "POS"]].astype(str)
    result_df['Site'] = vcf_file[['CHROM', 'POS']].agg(''.join, axis=1)
    columns_to_split = vcf_file.columns[9:]
    
    # This creates a for loop where we split the columns and store in a different data frame
    new_vcfs = []
    for col in columns_to_split:
        new_vcf = vcf_file[col].astype(str).str.split("|", expand=True)
        new_vcf.columns = [f"{col}_1", f"{col}_2"]
        new_vcfs.append(new_vcf)
    # Concatenate the new DataFrames along with the old data frame
    vcf_file = pd.concat([vcf_file] + new_vcfs, axis=1)
    # Drop the original columns that were split, leaving only split columns
    vcf_file.drop(columns=columns_to_split, inplace=True)
    # This calculates the allele frequencies of the parents
    allele_freq = np.empty([len(vcf_file), 2], dtype=float)
    for i in range(0,len(vcf_file)):
        if (i != 0):
            allele_freq[i-1,0] = round(parent_tot/4,2)
        parent_tot = 0
        for j in range(9,13):
            parent_tot = float(vcf_file.iloc[i,j]) + parent_tot
    allele_freq[i,0] = round(parent_tot/4,2)
    # This loop will calculate the parcentage in the progeny
    prog_num = len(vcf_file.columns) - 13
    # print(allele_freq)
    for i in range(0,len(vcf_file)):
            if (i != 0):
                allele_freq[i-1,1] = round(prog_tot,2)
            prog_tot = 0
            for j in range(14,len(vcf_file.columns)):
                prog_tot = float(vcf_file.iloc[i,j]) + prog_tot
    allele_freq[i,1] = round(prog_tot,2)
    result_df2 = pd.DataFrame(allele_freq, columns=['Parents', 'Progeny'])
    result_df = pd.concat([result_df, result_df2.reindex(result_df.index)], axis=1)
    print(result_df)
    result_df.to_csv(freq_file)
    print(prog_num)
    

    # parent files are at 9-12
    # 13: will be all other files

def VCF_importer(vcf):
    with open(vcf, 'r') as f:
        lines = [l for l in f if not l.startswith('##')]
    vcf_file = pd.read_csv(
        io.StringIO(''.join(lines)),
        dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
               'QUAL': str, 'FILTER': str, 'INFO': str},
        sep='\t'
    ).rename(columns={'#CHROM': 'CHROM'})
    return(vcf_file)

def csv_importer(file):
    data = pd.read_csv(file)
    return(data)



def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-vcf", "--vcf", type=str, help="A vcf file")
    parser.add_argument("-p", "--parent", type=str, help="A csv file containing parental data")
    parser.add_argument("--debug", default=False, action="store_true")
    args = parser.parse_args()

    # Handle debug flag
    global debug
    debug = args.debug
    return parser.parse_args()



main()