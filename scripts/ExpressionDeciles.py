import pandas as pd
import numpy as np

def convert_abundance(tsv_file, output_prefix):
    # Load the tx2gene file (transcript to gene mapping)
    tx2gene_file = "Homo_sapiens.GRCh38.Ensembl104.tx2gene"
    tx2gene = pd.read_csv(tx2gene_file, sep="\t", header=None, names=["transcript_id", "gene_id"])

    # Load transcript-level abundance file (from kallisto)
    abundance_data = pd.read_csv(tsv_file, sep="\t")
    
    # Merge with tx2gene to assign gene IDs to transcripts (using 'target_id' from abundance data)
    abundance_data = abundance_data.merge(tx2gene, how='left', left_on='target_id', right_on='transcript_id')

    # Summing abundance over genes, grouping by gene_id
    gene_abundance = abundance_data.groupby('gene_id')['tpm'].sum().reset_index()

    # Separate genes with zero abundance
    df0 = gene_abundance[gene_abundance['tpm'] == 0].copy()  # Ensure copy
    df1 = gene_abundance[gene_abundance['tpm'] != 0].copy()  # Ensure copy
    
    # Standardize non-zero abundance values (scale)
    standardized_col = (df1['tpm'] - df1['tpm'].mean()) / df1['tpm'].std()

    # Compute quantiles (deciles) for the standardized abundance values
    deciles = np.percentile(standardized_col, np.arange(0, 110, 10))

    # Assign deciles
    df1.loc[:, 'deciles'] = pd.cut(standardized_col, bins=deciles, labels=False, include_lowest=True) + 1
    
    # Set deciles for genes with zero abundance
    df0.loc[:, 'deciles'] = 1

    # Combine the zero-abundance and non-zero-abundance data
    output = pd.concat([df1[['gene_id', 'tpm', 'deciles']], df0[['gene_id', 'tpm', 'deciles']]])

    # Rename the 'tpm' column to 'abundance_TPM'
    output.rename(columns={'tpm': 'abundance_TPM'}, inplace=True)

    # Sort by gene_id
    output = output.sort_values(by='gene_id')

    # Ensure numeric columns are correctly formatted (avoid unwanted string conversion)
    output['abundance_TPM'] = pd.to_numeric(output['abundance_TPM'], errors='coerce')  # Force to float
    output['deciles'] = output['deciles'].astype(int)

    # Output file using pandas without unwanted formatting or quotes
    output.to_csv(f'{output_prefix}.expression_deciles_kallisto50_ensembl104.tsv', sep="\t", index=False, header=True, quoting=False)

if __name__ == "__main__":
    import sys
    tsv_file = sys.argv[1]
    output_prefix = sys.argv[2]
    convert_abundance(tsv_file, output_prefix)
