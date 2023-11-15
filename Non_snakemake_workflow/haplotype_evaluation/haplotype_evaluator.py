#! Python
# Objective is to load in a vcf file as a data frame and create tables that can load into R for visulization

# test file
# /home/drt06/Documents/Tall_fescue/Mapping_and_QTL/Mapping_and_QTL/Data/Real_Data/VCF/haplotypes/beagle_315x320.vcf
# /home/drt06/Documents/Tall_fescue/Mapping_and_QTL/Mapping_and_QTL/Data/Real_Data/Lists/usable_predicted_parents_double.csv
# /home/drt06/Documents/Tall_fescue/Mapping_and_QTL/Mapping_and_QTL/Data/Real_Data/VCF/haplotypes/Haplotype_Eval_Data

import pandas as pd
import argparse
import io
import os


def main():
    args = parse_args()
    vcf_file = VCF_importer(args.vcf)
    parent_file = csv_importer(args.parent)
    
    haplotype_percents(vcf_file, parent_file)



def haplotype_percents(vcf_file, parent_file):
    # Define window and step size
    window_size = 15
    step = 5

    # Define Variables
    data_folder = "/home/drt06/Documents/Tall_fescue/Mapping_and_QTL/Mapping_and_QTL/Data/Real_Data/VCF/haplotypes/Haplotype_Eval_Data"
    plant_names = vcf_file.columns.tolist()
    tot_prog = len(parent_file.index)

    # This loop iterates through the parent file to create variables for the current working individuals
    for prog_num in range (0, 3):
        current_proj = parent_file.iloc[prog_num]["Progeny"]
        current_parent1 = parent_file.iloc[prog_num]["Parent1"]
        current_parent2 = parent_file.iloc[prog_num]["Parent2"]
        print(current_proj,current_parent1, current_parent2 )
    
    test_prog = "315-4-43"
    test_parent1 = "320"
    test_parent2 = "315"
    # This loop checks if progeny profile is present in the vcf then starts to work with it and its parents.
    if test_prog in plant_names:
        print(f"{test_prog} is in the list.")
        result_df = pd.DataFrame(columns=['ProgAxParent1A', 'ProgAxParent1B', 'ProgAxParent2A', 'ProgAxParent2B', 'ProgBxParent1A', 'ProgBxParent1B', 'ProgBxParent2A', 'ProgBxParent2B'])
        percent_file = os.path.join(data_folder, test_prog)

        for i in range(0, len(vcf_file) - window_size + 1, step):
            # Extract windows from each column
            window_prog = vcf_file[test_prog].iloc[i:i+window_size]
            window_parent1 = vcf_file[test_parent1].iloc[i:i+window_size]
            window_parent2 = vcf_file[test_parent2].iloc[i:i+window_size]
            print(window_prog)
            # splitting windows
            prog_split_data = window_prog.str.split("|", expand=True)
            prog_split_data.columns = ['ProgA', 'ProgB']
            parent1_split_data = window_parent1.str.split("|", expand=True)
            parent1_split_data.columns = ['Parent1A', 'Parent1B']
            parent2_split_data = window_parent2.str.split("|", expand=True)
            parent2_split_data.columns = ['Parent2A', 'Parent2B']

            # Doing comparisons between sliding windows and saving the results to a data frame 
            percentage_ProgAxParent1A = (prog_split_data['ProgA'] == parent1_split_data['Parent1A']).mean()
            percentage_ProgAxParent1B = (prog_split_data['ProgA'] == parent1_split_data['Parent1B']).mean()
            percentage_ProgAxParent2A = (prog_split_data['ProgA'] == parent2_split_data['Parent2A']).mean()
            percentage_ProgAxParent2B = (prog_split_data['ProgA'] == parent2_split_data['Parent2B']).mean()
            percentage_ProgBxParent1A = (prog_split_data['ProgB'] == parent1_split_data['Parent1A']).mean()
            percentage_ProgBxParent1B = (prog_split_data['ProgB'] == parent1_split_data['Parent1B']).mean()
            percentage_ProgBxParent2A = (prog_split_data['ProgB'] == parent2_split_data['Parent2A']).mean()
            percentage_ProgBxParent2B = (prog_split_data['ProgB'] == parent2_split_data['Parent2B']).mean()
            
            new_data = pd.DataFrame({
                'ProgAxParent1A': [percentage_ProgAxParent1A],
                'ProgAxParent1B': [percentage_ProgAxParent1B],
                'ProgAxParent2A': [percentage_ProgAxParent2A],
                'ProgAxParent2B': [percentage_ProgAxParent2B],

                'ProgBxParent1A': [percentage_ProgBxParent1A],
                'ProgBxParent1B': [percentage_ProgBxParent1B],
                'ProgBxParent2A': [percentage_ProgBxParent2A],
                'ProgBxParent2B': [percentage_ProgBxParent2B]
            })
            # Concatenate the existing DataFrame with the new data
            result_df = pd.concat([result_df, new_data], ignore_index=True)

            result_df.to_csv(percent_file)

    else:
        print("progeny is missing") 
    


    # loop for nrows (maybe for now just read the ones on the first row)
    # read the progeny, read the parent1, read the parent2
    # create file called Progeny_name_hap1 and another file called Progeny_name_hap2
    # subset vcf file to only have progeny, parent1, and parent2
    # create a sliding window and save to a bunch of variables. RP1 RP2 Rp, LP1 LP2 Lp
    # Compare each Xp to all RP, whichever it is most similar to put that in file Progeny_name_hapx
    



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