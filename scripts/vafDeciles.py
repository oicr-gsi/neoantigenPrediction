import pandas as pd
import numpy as np
import sys

def calculate_vaf_deciles(variants_file, output_prefix):
    # Step 1: Extract the relevant fields from the VCF file using bcftools (similar to the original bash command)
    # Here, assume that `bcftools query` has already been run to generate the "ID.all_variants_vaf.tsv" file
    
    # Read in the VCF file (VAF values have been already extracted by `bcftools query` command)
    data = pd.read_csv(variants_file, sep="\t", header=None, names=["chr", "pos", "ref", "alt", "vaf"])

    # Step 2: Standardize the VAF column (z-score normalization)
    standardized_col = (data['vaf'] - data['vaf'].mean()) / data['vaf'].std()

    # Step 3: Calculate the deciles for the standardized VAF values
    deciles = standardized_col.quantile(np.arange(0, 1.1, 0.1)).values


    # Step 4: Assign deciles
    data['decile'] = pd.cut(standardized_col, bins=deciles, labels=False, include_lowest=True) + 1

    # Step 5: Save the results to a new file
    data.to_csv(f'{output_prefix}.deciles.tsv', sep="\t", index=False, header=True, quoting=False)


if __name__ == "__main__":
    # The VCF file is passed as a command line argument
    variants_file = sys.argv[1]
    output_prefix = sys.argv[2]
    calculate_vaf_deciles(variants_file, output_prefix)
