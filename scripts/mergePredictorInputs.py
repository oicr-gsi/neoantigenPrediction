import sys
import pandas as pd
import os

def merge_predictor_inputs(variants_peptides, variant_deciles, expression_deciles, rnaseq_variants):
    # Check if input files exist
    for file in [variants_peptides, variant_deciles, expression_deciles, rnaseq_variants]:
        if not os.path.exists(file):
            print(f"Error: {file} does not exist.")
            sys.exit(1)
    
    # Read in the TSV/CSV files
    vafPeptides = pd.read_csv(variants_peptides, sep="\t", header=0)
    vafDeciles = pd.read_csv(variant_deciles, sep="\t", header=0)
    vafDeciles.rename(columns={vafDeciles.columns[-1]:"vaf_decile"}, inplace=True)
    
    expressionDeciles = pd.read_csv(expression_deciles, sep="\t", header=0)
    expressionDeciles.rename(columns={expressionDeciles.columns[-1]:"expression_decile"}, inplace=True)
    
    variantsRNASeq = pd.read_csv(rnaseq_variants, sep="\t", header=None)
    variantsRNASeq.columns = ["chr", "pos", "ref", "alt", "variant_in_rna"]
    
    # Ensuring the 'pos' column is of type str in all DataFrames
    vafPeptides['pos'] = vafPeptides['pos'].astype(str)
    vafDeciles['pos'] = vafDeciles['pos'].astype(str)
    expressionDeciles['gene_id'] = expressionDeciles['gene_id'].astype(str)
    variantsRNASeq['pos'] = variantsRNASeq['pos'].astype(str)
    
    # Merge VAF Peptides with VAF Deciles
    m = pd.merge(vafPeptides, vafDeciles, on=["chr", "pos", "ref", "alt"], how="left")
    
    # Merge with Expression Deciles on Ensembl ID
    m = pd.merge(m, expressionDeciles, left_on="ensembl_id", right_on="gene_id", how="left")
    
    # Fill missing values in expression_decile with 1
    m["expression_decile"] = m["expression_decile"].fillna(1)
    
    # Merge with RNA-seq variants
    m = pd.merge(m, variantsRNASeq, on=["chr", "pos", "ref", "alt"], how="left")
    
    # Fill missing RNA-seq variant information with 0
    m["variant_in_rna"] = m["variant_in_rna"].fillna(0)
    
    
    
    # Defining the column order for the .tsv file
    columns = [
        'chr', 'pos', 'ref', 'alt', 'ensembl_id', 'rs_id', 'type', 'callers', 
        'tumor_dp', 'normal_dp', 'tumor_vaf', 'normal_vaf', 'consequence', 'impact', 
        'gene', 'transcript_id', 'biotype', 'nuc_change', 'aa_change', 'cdna_pos', 
        'cds_pos', 'protein_pos', 'amino_acids', 'codons', 'wt_peptides', 'mt_peptides', 
        'vaf', 'vaf_decile', 'abundance_TPM', 'expression_decile', 'variant_in_rna'
    ]
    
    # Making sure the columns in the DataFrame match the desired order
    m = m[columns]
    
    # Write the DataFrame to a TSV file
    m.to_csv('ID.output.merged.tsv', sep='\t', index=False, quoting=False)
    
    # Create a new DataFrame with formatted columns
    df = pd.DataFrame({
        'Unique identifier': m['chr'] + ";" + m['pos'] + ";" + m['ref'] + ";" + m['alt'],
        'Wt nmer': m['wt_peptides'],
        'Mut nmer': m['mt_peptides'],
        'Exome VAF decile': m['vaf_decile'],
        'Gene expression decile': m['expression_decile'],
        'Present in RNA-seq data': m['variant_in_rna']
    })

    # Write the new DataFrame to an Excel file
    df.to_excel('ID.output.merged.xls', index=False, engine='openpyxl')
    
    return 'ID.output.merged.tsv', 'ID.output.merged.xls'

if __name__ == "__main__":
    variants_peptides = sys.argv[1]
    variant_deciles = sys.argv[2]
    expression_deciles = sys.argv[3]
    rnaseq_variants = sys.argv[4]
    
    merge_predictor_inputs(variants_peptides, variant_deciles, expression_deciles, rnaseq_variants)
