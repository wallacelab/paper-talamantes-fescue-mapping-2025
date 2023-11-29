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
    window_size = 15
    step = 3
    # args.vcf is the vcf file
    # args.parent is a file of predicted parents
    args = parse_args()
    vcf_file = VCF_importer(args.vcf)
    parent_file = csv_importer(args.parent)
    # haplotype_percents(vcf_file, parent_file, data_folder, window_size, step)
    # haplotype_percents_parents(vcf_file, parent_file, data_folder, window_size, step)
    allelfrequncie(vcf_file,data_folder)





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