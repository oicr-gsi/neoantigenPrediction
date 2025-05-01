version 1.0

struct GenomeResources {
  String refModule
  String refFasta
}

struct VariantCalls {
  File vcf
  File vcfIndex
  String tumorId
}

struct HLACalls {
  Array[File] files
  Array[String] callers
}

workflow neoantigenPredictor {
  input {
    HLACalls HLAFiles
    VariantCalls DNAVariantCalls
    VariantCalls RNAVariantCalls
    File RNAAbundance
    String reference = "hg38"
    String outputFilePrefix
  }

  parameter_meta {
    HLAFiles : "an array of text files from all the HLA predictions with an array of the corresponding caller names"
    DNAVariantCalls : "the ensemble/combined DNA vcf file, from multiple callers, with the index and the tumourID"
    RNAVariantCalls : "the RNA seq variant calls from Haplotype Caller"
    RNAAbundance : "the expression data in text format"
    reference : "the reference build, defaults to hg38"
    outputFilePrefix : "a prefix to add to each of the output files, generally identifying the sample being analyzed"
  }

  meta {
    author : "Lawrence Heisler"
    email : "lheisler@oicr.on.ca"
    description : "A workflow that will use variant calls, expression data and HLA typing to predict Neoantigens"
    dependencies : [ 
      {
        name : "PCGR : Personal Cancer Genome Reporter, version 2.0.3",
        url : "https://github.com/sigven/pcgr/tree/v2.0.3"
      },
      {
        name : "bcftools version 1.9",
        url : "https://github.com/samtools/bcftools"
      },
      {
        name : "Variant Effect Predictor, version 112.0",
        url : "https://github.com/Ensembl/ensembl-vep"
      },
      {
        name : "pvactools version 4.3.0",
        url : "https://github.com/griffithlab/pVACtools"
      },
      {
        name : "SB_neoantigen_Models",
        url : "https://github.com/JaredJGartner/SB_neoantigen_Models"
      }
    ]
    output_meta: {
      NeoAntigenPredictions:{
        description: "An xlsx file with worksheets detailing the predictions",
        vidarr_label: "NeoAntigenPredictions"
      },
      NeoAntigenNmers:{
        description: "a tsv file with the predictions",
        vidarr_label: "NeoAntigenNmers"
      }
    }
  }

  Map[String,GenomeResources] resources = {
    "hg19": {
      "refModule": "hg19/p13", 
      "refFasta": "$HG19_ROOT/hg19_random.fa"
    },
    "hg38": {
      "refModule": "hg38/p12",
      "refFasta": "$HG38_ROOT/hg38_random.fa"
    }
  }
  
  output {
    File NeoAntigenPredictions = mergePredictorOutputs.xlsx
    File NeoAntigenNmers = mergePredictorOutputs.tsv
  }
  ### parse the HLA outputs to construct a string of HLAs
  call extractHLAs {
    input:
      hlafiles = HLAFiles.files,
      hlacallers = HLAFiles.callers,
      outputFilePrefix = outputFilePrefix
  }
  ### prepare for PCGR by adding in required INFO fields : TDP,TVAF,NDP,NVAF
  call format2pcgr {
    input:
      vcfin = DNAVariantCalls.vcf,
      tumorId = DNAVariantCalls.tumorId,
      outputFilePrefix = outputFilePrefix
  }
  ### call the PCGR software to determine which variants to keep as candidate sites
  call PCGR {
    input:
      vcf = format2pcgr.vcfout,
      vcfIndex = format2pcgr.vcfoutIndex,
      outputFilePrefix = outputFilePrefix
  }
 call vepAnnotate {
    input:
      vcf = PCGR.candidateCalls,
      refFasta = resources[reference].refFasta,
      outputFilePrefix = outputFilePrefix
  }
  call getPeptides {
    input: 
      vcf = vepAnnotate.annotatedCandidateCalls,
      tumorId = DNAVariantCalls.tumorId,
      outputFilePrefix = outputFilePrefix
  }
  call formatCalls {
    input:
      peptides = getPeptides.peptides,
      vcf = vepAnnotate.annotatedCandidateCalls,
      outputFilePrefix = outputFilePrefix
  }
  ### generate Deciles from the formatted vcf file
  call vafDeciles {
    input: 
      vcf= format2pcgr.vcfout,
      outputFilePrefix = outputFilePrefix
  }
  call ExpressionDeciles {
    input:
      tsv = RNAAbundance,
      outputFilePrefix = outputFilePrefix
  }
  call rnaseqVariants {
    input:
      vcf = RNAVariantCalls.vcf,
      outputFilePrefix = outputFilePrefix
  }
  call mergePredictorInputs {
     input:
       variants_peptides = formatCalls.tsv,
       variant_deciles = vafDeciles.deciles,
       expression_deciles = ExpressionDeciles.deciles,
       rnaseq_variants = rnaseqVariants.tsv,
       outputFilePrefix = outputFilePrefix
  }
  call chunkPredictorInputFile {
    input:
      xls = mergePredictorInputs.xls,
      chunksize = 30
  }
  scatter(predictorInput in chunkPredictorInputFile.predictorInputs){
    call predict {
       input:
         xls = predictorInput,
         hlas = extractHLAs.hlas
    }
  }
  call mergePredictorOutputs {
   input:
     predictorOutputs = predict.predictorOutput,
     predictorInputTSV = mergePredictorInputs.tsv,
     outputFilePrefix = outputFilePrefix
  }
}

task extractHLAs{
  input{
    Array[File] hlafiles
    Array[String] hlacallers
    String outputFilePrefix
    String modules = "neopipe/1.1.0" 
    Int jobMemory = 6
    Int timeout = 20	
  }
  
  parameter_meta {
      hlafiles: "A comma separated list of hla outputs from supported tools (t1k, optitype)"
	  hlacallers: "A comma separated list of hla caller names outputs from supported tools (t1k, optitype), isn the same order as hlafiles"
      outputFilePrefix: "The prefix to use for the output files"
	  modules: "Names and versions of modules"
	  jobMemory: "Memory allocated for task in GB"
      timeout: "Timeout in hours"
  }
  
  command<<<
  
  python3 $NEOPIPE_ROOT/bin/extractHLAs.py \
    "~{sep=" " hlafiles}" \
    "~{sep=" " hlacallers}" \
    ~{outputFilePrefix}.hlastring.txt

  >>>

  runtime {
    memory:  "~{jobMemory} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }
  
  output {
    String hlas = read_string("~{outputFilePrefix}.hlastring.txt")
  } 
}


task format2pcgr{
  input {
    File vcfin
    String tumorId
    String outputFilePrefix
    String modules = "neopipe/1.1.0 bcftools/1.9"
    Int jobMemory = 6
    Int timeout = 20
  }
  parameter_meta {
      vcfin: "dna variant calls in vcf format"
	  tumorId: "a string to identify the tumour column in the vcf file"
      outputFilePrefix: "The prefix to use for the output files"
	  modules: "Names and versions of modules"
	  jobMemory: "Memory allocated for task in GB"
      timeout: "Timeout in hours"
  }
  command <<<
  
  python3 $NEOPIPE_ROOT/bin/format2pcgr.py \
        -i ~{vcfin} \
        -o temp.ensemble.somatic.vt.annot.2callers.vcf.gz \
        -f 2 \
        -t ~{tumorId} \
        -v somatic \
        > format2pcgr.log 2>&1

  bcftools view -Oz -i 'TDP>=10 && TVAF>=0.05 && NDP>=10 && NVAF<=0.02' \
        -o ~{outputFilePrefix}.ensemble.somatic.vt.annot.2callers.vcf.gz \
        temp.ensemble.somatic.vt.annot.2callers.vcf.gz
  tabix -pvcf ~{outputFilePrefix}.ensemble.somatic.vt.annot.2callers.vcf.gz
  
  >>>
  runtime {
    memory:  "~{jobMemory} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }
  output {
    File vcfout = "~{outputFilePrefix}.ensemble.somatic.vt.annot.2callers.vcf.gz"
    File vcfoutIndex = "~{outputFilePrefix}.ensemble.somatic.vt.annot.2callers.vcf.gz.tbi"
  }
}


task PCGR{
  input {
    File vcf
    File vcfIndex
    String outputFilePrefix
    String modules = "neopipe/1.1.0 pcgr/2.0.3"
    Int jobMemory = 6
    Int timeout = 20
  }
  parameter_meta {
      vcf: "dna variant calls in PCGR-ready vcf format"
	  vcfIndex: "the vcf index file"
      outputFilePrefix: "The prefix to use for the output files"
	  modules: "Names and versions of modules"
	  jobMemory: "Memory allocated for task in GB"
      timeout: "Timeout in hours"
  }
  command <<<

  set -euxo pipefail
  
  mkdir pcgr
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
  
  R --vanilla<<RCODE
    library(pcgrr)
    yaml_fname<-"pcgr/~{outputFilePrefix}.pcgr.grch38.conf.yaml"
    my_log4r_layout <- function(level, ...) {
       paste0(format(Sys.time()), " - pcgr-report-generation - ",level, " - ", ..., "\n", collapse = "")
    }
   log4r_logger <- log4r::logger(threshold = "INFO", appenders = log4r::console_appender(my_log4r_layout))
   ## this gets passed on to all the log4r_* functions inside the pkg
   options("PCGRR_LOG4R_LOGGER" = log4r_logger)
  
   ## Generate report content
   pcg_report <- pcgrr::generate_report(yaml_fname = yaml_fname)
   pcgrr::write_report_tsv(report = pcg_report, output_type = 'snv_indel')
   pcgrr::write_report_excel(report = pcg_report)
  RCODE
  
  ### identify sites of interest (exonic, or non-exonic but actionable)
  python3 $NEOPIPE_ROOT/bin/filter_snv_indel.py ~{outputFilePrefix}
  
  ### now filter based on sites
  bcftools view -R ~{outputFilePrefix}.pcgr_filter_sites.txt -o ~{outputFilePrefix}.candidate_sites.vcf ~{vcf}

  	
  >>>
  runtime {
    memory:  "~{jobMemory} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }
  output {
    File candidateCalls = "~{outputFilePrefix}.candidate_sites.vcf"
  }
}

task vepAnnotate{
  input{
    File vcf
    String refFasta
    String outputFilePrefix
    String modules = "vep/112.0 pcgr/2.0.3 pvactools/4.3.0 hg38/p12"
    Int jobMemory = 6
    Int timeout = 20
   }

   parameter_meta {
      vcf: "dnaseq variant data in vcf format"
	  refFasta: "the genome reference in fasta format"
      outputFilePrefix: "The prefix to use for the output files"
	  modules: "Names and versions of modules"
	  jobMemory: "Memory allocated for task in GB"
      timeout: "Timeout in hours"
  }

   command<<<
   
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
   >>>

  runtime {
    memory:  "~{jobMemory} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }
  output {
    File annotatedCandidateCalls = "~{outputFilePrefix}.vep.vcf"
  }
}

task getPeptides {
  input{
    File vcf
    String tumorId
    String outputFilePrefix
    String modules = "pvactools/4.3.0"
    Int jobMemory = 6
    Int timeout = 20
   }

   parameter_meta {
      vcf: "candidate dna variants in vcf format"
	  tumorId: "column name for the tumour sample"
      outputFilePrefix: "The prefix to use for the output files"
	  modules: "Names and versions of modules"
	  jobMemory: "Memory allocated for task in GB"
      timeout: "Timeout in hours"
   }  
   
   command<<<

   pvacseq generate_protein_fasta \
   -s ~{tumorId} \
   -d 12 \
   ~{vcf} 12 ~{outputFilePrefix}.peptides.fa
   
   ## convert to single line per site with WTid,WTseq,MTid,MTseq
   cat ~{outputFilePrefix}.peptides.fa | paste - - - - > ~{outputFilePrefix}.peptides.0.txt
   ## keep only records where MT != WT peptide
   cat ~{outputFilePrefix}.peptides.0.txt | awk '{if($2 != $4) print}' > ~{outputFilePrefix}.peptides.1.txt
   ## order by numeric identifier
   cat ~{outputFilePrefix}.peptides.1.txt | sort -t . -k2 -g > ~{outputFilePrefix}.peptides.2.txt
   ## keep only if unique and peptide length == 25
   cat ~{outputFilePrefix}.peptides.2.txt | awk -F'\t' '!seen[$2,$4]++ && length($2) >= 25 && length($4) >= 25' > ~{outputFilePrefix}.peptides.txt


   >>>

  runtime {
    memory:  "~{jobMemory} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }
  output {
    File peptides = "~{outputFilePrefix}.peptides.txt"
  }
}

task formatCalls{
  input {
    File peptides
    File vcf
    String outputFilePrefix
    String modules = "neopipe/1.1.0 bcftools/1.9"
    Int jobMemory = 6
    Int timeout = 20
  }

  parameter_meta {
      peptides: "peptides generated from the candidate calls"
	  vcf: "candidate calls in vcf format"
      outputFilePrefix: "The prefix to use for the output files"
	  modules: "Names and versions of modules"
	  jobMemory: "Memory allocated for task in GB"
      timeout: "Timeout in hours"
  }  

  command<<<

  python3 $NEOPIPE_ROOT/bin/formatCalls.py \
    ~{peptides} \
    ~{vcf} \
    ~{outputFilePrefix}
    
    >>>  

  runtime {
    memory:  "~{jobMemory} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }
  output {
    File tsv = "~{outputFilePrefix}.candidateCalls.tsv"
  }
}

task vafDeciles{
  input {
    File vcf
    String outputFilePrefix
    String modules = "neopipe/1.1.0 bcftools/1.9"
    Int jobMemory = 6
    Int timeout = 20
  }

  parameter_meta {
      vcf: "dnaseq variant data in vcf format"
      outputFilePrefix: "The prefix to use for the output files"
	  modules: "Names and versions of modules"
	  jobMemory: "Memory allocated for task in GB"
      timeout: "Timeout in hours"
  }
    
  command<<<

  ## extract relevant fields from the vcf records
  bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%TVAF\n' ~{vcf} > ~{outputFilePrefix}.all_variants_vaf.tsv
  
  ## generated deciles
  python3 $NEOPIPE_ROOT/bin/vafDeciles.py \
  ~{outputFilePrefix}.all_variants_vaf.tsv \
  ~{outputFilePrefix} \

  >>>  

  runtime {
    memory:  "~{jobMemory} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }
  output {
    File deciles = "~{outputFilePrefix}.deciles.tsv"
  }
}

task ExpressionDeciles{
  input {
    File tsv
    String outputFilePrefix
    String modules = "neopipe/1.1.0 ensembl/104-hg38"
    Int jobMemory = 6
    Int timeout = 20
  }

  parameter_meta {
      tsv: "A tsv file with expression data"
      outputFilePrefix: "The prefix to use for the output files"
	  modules: "Names and versions of modules"
	  jobMemory: "Memory allocated for task in GB"
      timeout: "Timeout in hours"
  }
  
  command<<<

  python3 $NEOPIPE_ROOT/bin/ExpressionDeciles.py \
    ~{tsv} \
    ~{outputFilePrefix}

  >>>  

  runtime {
    memory:  "~{jobMemory} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }
  output {
    File deciles = "~{outputFilePrefix}.expression_deciles_kallisto50_ensembl104.tsv"
  }
}


task rnaseqVariants{
  input {
    File vcf
    String outputFilePrefix
    String modules = "bcftools/1.9"
    Int jobMemory = 6
    Int timeout = 20
  }
  
  parameter_meta {
      vcf: "rnaseq variant data in vcf format"
      outputFilePrefix: "The prefix to use for the output files"
	  modules: "Names and versions of modules"
	  jobMemory: "Memory allocated for task in GB"
      timeout: "Timeout in hours"
  }
  command<<<
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t1\n' ~{vcf} > ~{outputFilePrefix}.rnaseq_variants.tsv
    >>>  

  runtime {
    memory:  "~{jobMemory} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }
  output {
    File tsv = "~{outputFilePrefix}.rnaseq_variants.tsv"
  }
}


task mergePredictorInputs{
   input {
     File variants_peptides
     File variant_deciles
     File expression_deciles
     File rnaseq_variants
     String outputFilePrefix
     String modules = "neopipe/1.1.0 sb-neoantigen-models/1.0.0"
     Int jobMemory = 6
     Int timeout = 20	  
   }

   parameter_meta {
      variants_peptides: "tsv file with candidate variants and associated peptide predictions"
	  variant_deciles: "vaf deciles for candidate variants"
      expression_deciles: "expression deciles for expressed sequences "
	  rnaseq_variants: "rnaseq variants"
	  modules: "Names and versions of modules"
	  jobMemory: "Memory allocated for task in GB"
      timeout: "Timeout in hours"
   }

   command<<<

   python3 $NEOPIPE_ROOT/bin/mergePredictorInputs.py \
    ~{variants_peptides} \
    ~{variant_deciles} \
    ~{expression_deciles} \
    ~{rnaseq_variants} \
    ~{outputFilePrefix}

   >>>  

  runtime {
    memory:  "~{jobMemory} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }
  output {
    File tsv = "~{outputFilePrefix}.output.merged.tsv"
    File xls = "~{outputFilePrefix}.output.merged.xls"
  }

}


task predict{
  input{
    File xls
    String hlas
    String modules = "sb-neoantigen-models/1.0.0" 
    Int jobMemory = 6
    Int timeout = 20	  	
  }

  parameter_meta {
      xls: "xls file with inputs to sb_neoantigen_predictor"
	  hlas: "an string with a list of hlas"
	  modules: "Names and versions of modules"
	  jobMemory: "Memory allocated for task in GB"
      timeout: "Timeout in hours"
  }

  command<<<
  ### set up split to parallelize job
  python $SB_NEOANTIGEN_MODELS_ROOT/src/GenerateScores.py ~{xls} ~{hlas}
  ### set up join to join results from parallelized jobs
  >>>
  
  String prefix = basename(xls)
  
  runtime {
    memory:  "~{jobMemory} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }
  output {
    File predictorOutput = "~{prefix}_scored.xlsx"
  }  
}


task chunkPredictorInputFile {
   input{
     File xls
     Int chunksize
     String modules = "neopipe/1.1.0 sb-neoantigen-models/1.0.0"
     Int jobMemory = 6
     Int timeout = 20 
   }
  parameter_meta {
      xls: "xls file with inputs to sb_neoantigen_predictor"
    chunksize: "the number of records to include in each chunk"
    modules: "Names and versions of modules"
    jobMemory: "Memory allocated for task in GB"
      timeout: "Timeout in hours"
  }

  command<<<
  
  python3 $NEOPIPE_ROOT/bin/chunkPredictorInputFile.py \
    ~{xls} \
    ~{chunksize} \

  >>>

  runtime {
    memory:  "~{jobMemory} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }

  output {
    Array[File] predictorInputs = glob("predict_input*.xls")
  } 
}  


task mergePredictorOutputs {
   input{
     Array[File] predictorOutputs
     File predictorInputTSV
     String outputFilePrefix
     String modules = "neopipe/1.1.0 sb-neoantigen-models/1.0.0"
     Int jobMemory = 6
     Int timeout = 20 
   }
   
  parameter_meta {
    predictorOutputs: "outputs from sb_neoantigen_predictor that need to be merged"
    predictorInputTSV: "the full input to the neoantigen_predictor"
    modules: "Names and versions of modules"
    jobMemory: "Memory allocated for task in GB"
    timeout: "Timeout in hours"
  }
  
  
  command<<<

  python3 $NEOPIPE_ROOT/bin/mergePredictorOutputs.py \
  "~{sep=" " predictorOutputs}" \
  ~{predictorInputTSV} \
  ~{outputFilePrefix}

  >>>   


  runtime {
    memory:  "~{jobMemory} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }
   
   output {
     File xlsx = "~{outputFilePrefix}_neoantigenPredictions.xlsx"
     File tsv =  "~{outputFilePrefix}_neoantigenPredictions.tsv"
   }
}
