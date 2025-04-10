import sys
import pandas as pd
import numpy as np

def main(predictor_outputs, predictor_input_tsv, outputFilePrefix):
    # Read the Excel files and concatenate the sheets for "InputData", "NmersScored", and "MmpsScored"
    data = {}
    for ds in ["InputData", "NmersScored", "MmpsScored"]:
        for predictor_output in predictor_outputs:
            df = pd.read_excel(predictor_output, sheet_name=ds)
            if ds not in data:
                data[ds] = df
            else:
                data[ds] = pd.concat([data[ds], df], axis=0)

    # Write the data to an Excel file
    writer = pd.ExcelWriter(f'{outputFilePrefix}_neoantigenPredictions.xlsx', engine='xlsxwriter')
    data["InputData"].to_excel(writer, sheet_name='InputData', index=False)
    data["NmersScored"].to_excel(writer, sheet_name='NmersScored', index=False)
    data["MmpsScored"].to_excel(writer, sheet_name='MmpsScored', index=False)
    writer.close()
    
    # Annotate the predictorInputTSV with information from the Excel file
    # NmersScored has 1 record per Unique identifier
    # MmpsScores had multiple records per Unique identifer, across a range of MMP model score
    # merge the two, matching each of the NmersScored mmp scores to the MMP model score, concatenating values if there are multiple matches
    df_nmers = data["NmersScored"]
    df_mmps = data["MmpsScored"]

    # Prepare to capture the list of lists for matching MMP scores
    list_mmps_match = []

    # to capture the list of lists, for conversion to a dataframe
    for index, row in df_nmers.iterrows():
        match = []
        uid = row['Unique identifier']
        match.append(uid)

        # Adjust precision for MMP scores
        score1 = round(row["mmp score 1"], 7)
        score2 = round(row["mmp score 2"], 7)

        for score in [score1, score2]:
            # Match on uid and score
            mmps_rows = df_mmps[(df_mmps['Unique identifier'] == uid) & (round(df_mmps['MMP model score'], 7) == score)]
            
            # Concatenate values if multiple rows are found
            mutant_peptides = ";".join(mmps_rows['Mutant peptide'])
            match.append(mutant_peptides)
            hlas = ";".join(mmps_rows['HLA'])
            match.append(hlas)

        # Add the matched result to the list
        list_mmps_match.append(match)

    # Convert the list to a dataframe
    df_mmps_match = pd.DataFrame(list_mmps_match, columns=["Unique identifier", "mmpScore1Pep", "mmpScore1HLA", "mmpScore2Pep", "mmpScore2HLA"])

    # Merge the new columns to the NmersScored dataframe
    df_nmers_extended = pd.merge(df_nmers, df_mmps_match, on="Unique identifier")

    # Split the Unique identifier into separate columns and adjust dtype for pos
    df_nmers_extended[["chr", "pos", "ref", "alt"]] = df_nmers_extended['Unique identifier'].str.split(";", expand=True)
    df_nmers_extended = df_nmers_extended.drop(["Unique identifier"], axis=1)
    df_nmers_extended["pos"] = df_nmers_extended["pos"].astype('int64')

    # Load the TSV file
    df_tsv = pd.read_csv(predictor_input_tsv, sep="\t")

    # Merge the Nmer information into the TSV
    df_tsv_extended = pd.merge(df_tsv, df_nmers_extended, on=["chr", "pos", "ref", "alt"])

    # Save the extended TSV to a file
    df_tsv_extended.to_csv(f'{output_file_prefix}_neoantigenPredictions.tsv', index=False, sep="\t")

if __name__ == "__main__":
    # Expecting arguments from the WDL file execution (passed as sys.argv)
    predictor_outputs = sys.argv[1].split()  # Space-separated list of Excel file paths
    predictor_input_tsv = sys.argv[2]  # predictor input TSV file
    output_file_prefix = sys.argv[3]  # output file prefix

    main(predictor_outputs, predictor_input_tsv, output_file_prefix)
