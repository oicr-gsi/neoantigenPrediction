# neoantigenPredictor

A workflow that will use variant calls, expression data and HLA typing to predict Neoantigens

## Overview

## Dependencies

* [PCGR : Personal Cancer Genome Reporter, version 2.0.3](https://github.com/sigven/pcgr/tree/v2.0.3)
* [bcftools version 1.9](https://github.com/samtools/bcftools)
* [Variant Effect Predictor, version 112.0](https://github.com/Ensembl/ensembl-vep)
* [pvactools version 4.3.0](https://github.com/griffithlab/pVACtools)
* [SB_neoantigen_Models](https://github.com/JaredJGartner/SB_neoantigen_Models)


## Usage

### Cromwell
```
java -jar cromwell.jar run neoantigenPredictor.wdl --inputs inputs.json
```

### Inputs

#### Required workflow parameters:
Parameter|Value|Description
---|---|---
`HLAFiles`|HLACalls|an array of text files from all the HLA predictions with an array of the corresponding caller names
`DNAVariantCalls`|VariantCalls|the ensemble/combined DNA vcf file, from multiple callers, with the index and the tumourID
`RNAVariantCalls`|VariantCalls|the RNA seq variant calls from Haplotype Caller
`RNAAbundance`|File|the expression data in text format
`outputFilePrefix`|String|a prefix to add to each of the output files, generally identifying the sample being analyzed


#### Optional workflow parameters:
Parameter|Value|Default|Description
---|---|---|---
`reference`|String|"hg38"|the reference build, defaults to hg38


#### Optional task parameters:
Parameter|Value|Default|Description
---|---|---|---
`extractHLAs.modules`|String|"neopipe/1.0.0"|Names and versions of modules
`extractHLAs.jobMemory`|Int|6|Memory allocated for task in GB
`extractHLAs.timeout`|Int|20|Timeout in hours
`format2pcgr.modules`|String|"neopipe/1.0.0 bcftools/1.9"|Names and versions of modules
`format2pcgr.jobMemory`|Int|6|Memory allocated for task in GB
`format2pcgr.timeout`|Int|20|Timeout in hours
`PCGR.modules`|String|"pcgr/2.0.3"|Names and versions of modules
`PCGR.jobMemory`|Int|6|Memory allocated for task in GB
`PCGR.timeout`|Int|20|Timeout in hours
`vepAnnotate.modules`|String|"vep/112.0 pcgr/2.0.3 pvactools/4.3.0 hg38/p12"|Names and versions of modules
`vepAnnotate.jobMemory`|Int|6|Memory allocated for task in GB
`vepAnnotate.timeout`|Int|20|Timeout in hours
`getPeptides.modules`|String|"pvactools/4.3.0"|Names and versions of modules
`getPeptides.jobMemory`|Int|6|Memory allocated for task in GB
`getPeptides.timeout`|Int|20|Timeout in hours
`formatCalls.modules`|String|"bcftools/1.9"|Names and versions of modules
`formatCalls.jobMemory`|Int|6|Memory allocated for task in GB
`formatCalls.timeout`|Int|20|Timeout in hours
`vafDeciles.modules`|String|"neopipe/1.0.0 bcftools/1.9"|Names and versions of modules
`vafDeciles.jobMemory`|Int|6|Memory allocated for task in GB
`vafDeciles.timeout`|Int|20|Timeout in hours
`ExpressionDeciles.modules`|String|"neopipe/1.0.0 ensembl/104-hg38"|Names and versions of modules
`ExpressionDeciles.jobMemory`|Int|6|Memory allocated for task in GB
`ExpressionDeciles.timeout`|Int|20|Timeout in hours
`rnaseqVariants.modules`|String|"bcftools/1.9"|Names and versions of modules
`rnaseqVariants.jobMemory`|Int|6|Memory allocated for task in GB
`rnaseqVariants.timeout`|Int|20|Timeout in hours
`mergePredictorInputs.modules`|String|"neopipe/1.0.0"|Names and versions of modules
`mergePredictorInputs.jobMemory`|Int|6|Memory allocated for task in GB
`mergePredictorInputs.timeout`|Int|20|Timeout in hours
`chunkPredictorInputFile.modules`|String|"neopipe/1.0.0 sb-neoantigen-models/1.0.0"|Names and versions of modules
`chunkPredictorInputFile.jobMemory`|Int|6|Memory allocated for task in GB
`chunkPredictorInputFile.timeout`|Int|20|Timeout in hours
`predict.modules`|String|"sb-neoantigen-models/1.0.0"|Names and versions of modules
`predict.jobMemory`|Int|6|Memory allocated for task in GB
`predict.timeout`|Int|20|Timeout in hours
`mergePredictorOutputs.modules`|String|"neopipe/1.0.0 sb-neoantigen-models/1.0.0"|Names and versions of modules
`mergePredictorOutputs.jobMemory`|Int|6|Memory allocated for task in GB
`mergePredictorOutputs.timeout`|Int|20|Timeout in hours


### Outputs

Output | Type | Description | Labels
---|---|---|---
`NeoAntigenPredictions`|File|An xlsx file with worksheets detailing the predictions|vidarr_label: NeoAntigenPredictions
`NeoAntigenNmers`|File|a tsv file with the predictions|vidarr_label: NeoAntigenNmers


## Commands
 
 This section lists command(s) run by neoantigenPrediction workflow
 
 * Extraction of HLA Calls to an HLA String as input to sb_neoantigen_prediction
 
 ```
 	#Inline python code, see WDL task for extractHLAs for details (todo: convert to script)
 
 ```
 
 * Format the Variant Call input so that they are PCGR-ready.  This will keep only calls identified by 2 or more callers and filters on Depth and VAF
 
 ```
   # Run format2pcgr
   format2pcgr \
         -i ~{vcfin} \
         -o temp.ensemble.somatic.vt.annot.2callers.vcf.gz \
         -f 2 \
         -t ~{tumorId} \
         -v somatic \
         > format2pcgr.log 2>&1
 
   # Filter with vcftools
   bcftools view -Oz -i 'TDP>=10 && TVAF>=0.05 && NDP>=10 && NVAF<=0.02' \
         -o ~{outputFilePrefix}.ensemble.somatic.vt.annot.2callers.vcf.gz \
         temp.ensemble.somatic.vt.annot.2callers.vcf.gz
   
   # Index the vcf output
   tabix -pvcf ~{outputFilePrefix}.ensemble.somatic.vt.annot.2callers.vcf.gz
 
 ```
 
 
 * Run PCGR (Personal Cancer Genome Reporter)
 
 ```
    # Extract Calls + Tumour VAF with bcftools
    set -euxo pipefail
   
    mkdir pcgr
    ### pcgr is called without reporting, to generate data files from which annotated candidate sites can be extracted
    pcgr --vcf2maf \
     --force_overwrite \
     --vep_buffer_size 500 \
     --vep_regulatory \
     --exclude_dbsnp_nonsomatic \
     --exclude_likely_het_germline \
     --exclude_likely_hom_germline \
     --exclude_nonexonic \
     --maf_gnomad_global 0.002 \
     --tumor_site 0 \
     --assay WGS \
     --tumor_dp_tag TDP --tumor_af_tag TVAF --tumor_dp_min 10 --tumor_af_min 0.05 \
     --control_dp_tag NDP --control_af_tag NVAF --control_dp_min 10 --control_af_max 0.02 \
     --input_vcf ~{vcf} \
     --output_dir pcgr \
     --genome_assembly grch38 \
     --sample_id ~{outputFilePrefix} \
     --vep_dir $VEP_DIR \
     --refdata_dir $REFDATA_DIR \
     --pcgrr_conda $PCGR_ROOT \
     --no_reporting
 
    # Inline R code run pcgr reporting functions to generated tsv and excel files, see WDL task PCGR for details (todo: convert to script)
 
    # Inline python code that filters PCGR tsv output based on EXONIC_STATUS and ACTIONABILITY_TIER to generate candidate sites, see WDL task PCGR for details (todo: convert to script)
 
    # filter the input vcf to only the candidate sites
    bcftools view -R ~{outputFilePrefix}.pcgr_filter_sites.txt -o ~{outputFilePrefix}.candidate_sites.vcf ~{vcf}
 
 ```
 
 * Annotate the Candidate sites with VEP to add in consequence information that is needed for Peptide identification
 
 ```
 
    # Extract Calls + Tumour VAF with bcftools
    vep \
      --input_file ~{vcf} \
      --output_file ~{outputFilePrefix}.vep.vcf \
      --format vcf --vcf --symbol --terms SO --tsl \
      --hgvs --fasta ~{refFasta} \
      --offline \
      --cache --dir_cache $VEP_DIR --cache_version 112 \
      --plugin Frameshift --plugin Wildtype \
      --dir_plugins $PVACTOOLS_VEP_PLUGINS --pick --transcript_version \
      --force_overwrite
 
 ```
 
 
 * Generate Peptides for the Candidate sites with pVactools
 
 ```
    pvacseq generate_protein_fasta \
      -s ~{tumorId} \
      -d 12 \
      ~{vcf} 12 ~{outputFilePrefix}.peptides.fa
 
 ```
 
 
 * Combined candidate calls with peptide predictions for input to sb_neoantigen_prediction
 
 ```
    # Inline Python code, see WDL task formatCalls for details (todo: convert to script)
 
 ```
 
 
 * Generate Deciles from the PCGR-ready variant calls
 
 ```
 
    # Inline R code to generate Deciles, see WDL task ExpressionDeciles for details (todo: convert to script)
 
 ```
 
 * Generate Deciles from the rnaseq expression data
 
 ```
    # Extract Calls + Tumour VAF with bcftools
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%TVAF\n' ~{vcf} > ~{outputFilePrefix}.all_variants_vaf.tsv
 
    # Inline R code to generate Deciles, see WDL task vafDeciles for details (todo: convert to script)
 
 ```
 
 
 * Extract RNASeq variants to a tsv format
 
 ```
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t1\n' ~{vcf} > ~{outputFilePrefix}.rnaseq_variants.tsv
 
 ```
 
 
 * Combine all sb_neoantigen predictor inputs (candidate calls, vaf deciles, rnaseq deciles, expression calls) to a tsv and xlsx format
 
 ```
    # Inline R code to combine data, see WDL task mergePredictoriInputs for details (todo: convert to script)
 
 ```
 
 
 
 * process data through sb_neoantigen_predictor
 
 ```
    # Input is the xls file with combined input and the hla string
    python $SB_NEOANTIGEN_MODELS_ROOT/src/GenerateScores.py ~{xls} ~{hlas}
 
 ```
 
 * scatter/gather for predictor.  For improved performance the candidate inputs are separated into smaller subsets and sb_neoantigen_predictor is run in parallel.  The output from each subset is then merged into the final output
 
 ```
    # subsetting of the data is done with inline python code, see the WDL task chunkPredictorInputFile for details (todo: convert to script)
 
    # merging the subset outputs is done with inline python code, see the WDL task mergePredictorOutputs for details (todo: convert to script)
 
  
 
 ```
 
 
 
 
 
 ## Support

For support, please file an issue on the [Github project](https://github.com/oicr-gsi) or send an email to gsi@oicr.on.ca .

_Generated with generate-markdown-readme (https://github.com/oicr-gsi/gsi-wdl-tools/)_
