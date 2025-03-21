import re
import os
import sys

def main(peptides_file, vcf_file, output_prefix):
    # Ensure peptides_file is a string path
    if not isinstance(peptides_file, str):
        raise TypeError(f"Expected string for peptides_file, got {type(peptides_file)}")

    fout = open(f"{output_prefix}.candidateCalls.tsv", "w")
    
    output_fields = [
        "chr", "pos", "rs_id", "ref", "alt", "type", "callers", "tumor_dp", 
        "normal_dp", "tumor_vaf", "normal_vaf", "consequence", "impact", 
        "gene", "ensembl_id", "transcript_id", "biotype", "nuc_change", 
        "aa_change", "cdna_pos", "cds_pos", "protein_pos", "amino_acids", 
        "codons", "wt_peptides", "mt_peptides"
    ]
    print(*output_fields, sep="\t", file=fout)
    
    # Load and catalogue the peptides
    peptides = {}
    with open(peptides_file) as pfile:
        for line in pfile:
            line = line.rstrip("\n")
            WTid, WT_peptide, MTid, MT_peptide = line.split("\t")
            d = dict(zip(["prefix", "ordinal", "gene", "tid0", "tid1", "type", "change"], WTid.split(".")))
            d['wt_peptide'] = WT_peptide
            d['mt_peptide'] = MT_peptide
            d['transcript_id'] = d['tid0'] + "." + d['tid1']
            results = re.search("(\d+)(.*)", d['change'])
            d['amino_acids'] = results[2]
            d['protein_pos'] = results[1]
            key = "|".join([d[x] for x in ['gene', 'transcript_id', 'protein_pos', 'amino_acids']])
            peptides[key] = d
    
    # Keys/fields of interest
    info_keys = ["CHROM", "POS", "ID", "REF", "ALT", "TYPE", "CALLERS", "TDP", "NDP", "TVAF", "NVAF"]
    csq_keys = ["Consequence", "IMPACT", "SYMBOL", "Gene", "Feature", "BIOTYPE", "HGVSc", "HGVSp", 
                "cDNA_position", "CDS_position", "Protein_position", "Amino_acids", "Codons", 
                "wt_peptide", "mt_peptide"]
    
    # Get CSQ fields from the header
    headers = os.popen(f'bcftools view -h "{vcf_file}"').readlines()
    CSQ_fields = []
    for header in headers:
        if "##INFO" in header:
            if "ID=CSQ" in header:
                CSQ_fields_string = re.search("Format: (.*)\"", header)
                CSQ_fields = CSQ_fields_string[1].split("|")
    
    # Get required fields from each record
    recs = os.popen(f'bcftools query -f "%CHROM\t%POS\t%ID\t%REF\t%ALT\t%TYPE\t%CALLERS\t%TDP\t%NDP\t%TVAF\t%NVAF\t%CSQ\n" "{vcf_file}"').readlines()
    
    for rec in recs:
        rec = rec.rstrip("\n")    
        # Get the values for this record and store in info_dict
        v = rec.split("\t")
        info_dict = dict(zip(info_keys + ["CSQ"], v))
        
        # Get the values for the CSQ part of the record, store in csq_dict
        v = info_dict['CSQ'].split("|")
        csq_dict = dict(zip(CSQ_fields, v))
        
        # For a key for the peptide dictionary
        peptide_key = "|".join([csq_dict[x] for x in ['SYMBOL', 'Feature', 'Protein_position', 'Amino_acids']])
        
        # Print out the record only if the key is in the peptide dictionary 
        if peptide_key in peptides:
            csq_dict['wt_peptide'] = peptides[peptide_key]['wt_peptide']
            csq_dict['mt_peptide'] = peptides[peptide_key]['mt_peptide']
            # Get the required info and csq values (including the peptides), and print to the output file
            info_values = [info_dict[x] for x in info_keys]
            csq_values = [csq_dict[x] for x in csq_keys] 
            all_values = info_values + csq_values
            print(*all_values, sep="\t", file=fout)
    
    fout.close()

if __name__ == "__main__":
    peptides_file = sys.argv[1]
    vcf_file = sys.argv[2]
    output_prefix = sys.argv[3]
    main(peptides_file, vcf_file, output_prefix)
