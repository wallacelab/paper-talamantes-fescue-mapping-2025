# Purpose: This code is to be used after Making_all_Residual_data_subsets.R. 
# It will take all the outputs from that and calculate the differences in the H value
# Then create a file and output the name of the file along with the first time the H value difference 
# was less than .05 and did not start at 0.

# Author: Darrian Talamantes

import os
import pandas as pd

def process_files(directory, output_file):
    results = []

    # Iterate over all files in the directory
    for filename in os.listdir(directory):
        if "data_1" in filename and filename.endswith(".csv"):
            filepath = os.path.join(directory, filename)
            
            # Read the CSV file
            df = pd.read_csv(filepath)
            
            # Use only the first half of the rows
            half_index = len(df) // 2
            df_half = df.iloc[:half_index]

            first_match_found = False
            last_row = None

            # Process differences in the "H" column
            for i in range(1, len(df_half) - 1):  # Loop through rows starting from the second
                diff = df_half.iloc[i + 1]["H"] - df_half.iloc[i]["H"]
                last_row = df_half.iloc[i + 1]  # Keep track of the last row

                if not first_match_found and diff != 0 and abs(diff) <= 0.05 and df_half.iloc[i + 1]["H"] != 0 and  df_half.iloc[i]["H"] != 0:
                    # Save the first match and break
                    results.append({
                        "File": filename.replace("_data_1.csv", ""),
                        "N": df_half.iloc[i]["N"],
                        "H": df_half.iloc[i]["H"]
                    })
                    first_match_found = True
                    break
            
            # If no match was found, add the last row
            if not first_match_found and last_row is not None:
                results.append({
                    "File": filename.replace("_data_1.csv", ""),
                    "N": last_row["N"],
                    "H": last_row["H"]
                })

    # Save results to a new CSV file
    results_df = pd.DataFrame(results)
    results_df.to_csv(output_file, index=False)
    print(f"Results saved to {output_file}")

# Example usage
directory = "/home/darrian/Desktop/UGA/Wallace_Lab/Mapping_and_QTL/Data/Heritability_Outputs"  # Put Directory path here
output_file = "/home/darrian/Desktop/UGA/Wallace_Lab/Mapping_and_QTL/Data/Heritability_Outputs/Final_Heritabilities.csv"
process_files(directory, output_file)