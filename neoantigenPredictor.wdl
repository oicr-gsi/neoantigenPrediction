version 1.0

struct GenomeResources {
    String refModule
    String refFasta
}

struct VariantCalls {
  File vcf
  File vcfIndex
}

workflow neoantigenPredictor {
  input {
    Array[File] HLACalls
    VariantCalls DNAVariantCalls
    VariantCalls RNAVariantCalls
    File Expression
    String reference = "hg38"
  }

  parameter_meta {
    HLACalls : "an array of text files from all the HLA predictions"
    VariantsCallSet : "the set of vcf files from each somatic variant callers, with the index and caller identified"
    Expression : "the expression data in text format"
    VariantsRNA : "the vcf files from the rna variant caller, with the index and caller identified"
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
    File NeoAntigenPredictions = format4PCGR.vcfout
  }

  meta {
    author : "Lawrence Heisler"
    email : "lheisler@oicr.on.ca"
    description : "A workflow that will use variant calls, expression data and HLA typing to predict Neoantigens"
    dependencies : [ 
      {
        name : "xxx",
        url : "xxx"
      }
    ]
    output_meta: {
      NeoAntigenPredictions:{
        description: "xxx",
        vidarr_label: "xxx"
      }
    }
  }

  
  call format4PCGR {
    input:
	   vcfin = DNAVariantCalls["vcf"]
  }
  call PCGR {
    input:
      vcf = format4PCGR.vcfout,
      vcfIndex = format4PCGR.vcfoutIndex	   
  }

  call vafDeciles {
    input: vcf= DNAVariantCalls["vcf"]
  }

#  call VEPAnnotated {
#    input:
#  }

#  call pVAC {
#    input:
#  }

#  call VariantsWithPeptides {
#    input:
#  }

#  call ExpressionDeciles {
#    input:
#
#  }
 
#  call Predict {
#    input:
#  }



}


task vafDeciles{
  input {
    File vcf
    String modules = "bcftools/1.9"
    Int jobMemory = 6
    Int timeout = 20
  }
  
  command<<<
  module load bcftools
  ## extract relevant fields from the vcf records
  bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%TVAF\n' ~{vcf} > ID.all_variants_vaf.tsv
  
  ## generated deciles with R code
  R --vanilla <<< RCODE
  library(dplyr)

  # Function to calculate deciles and output results
  calculate_deciles <- function(file_path) {
    # Read the file
    data <- read.delim(file_path, header = FALSE, sep = "\t", col.names = c("Column1", "Column2", "Column3", "Column4", "Column5"))
  
    # Standardize the 5th column
    standardized_col <- scale(data$Column5)
  
    # Calculate deciles
    deciles <- quantile(standardized_col, probs = seq(0, 1, 0.1))
    #print(head(data))
    # Output results
    output <- data.frame(data[, 1:5], Deciles = cut(standardized_col, breaks = deciles, labels = FALSE, include.lowest = TRUE))
    return(output)
  }
  
  # Call the function on the selected file
  output <- calculate_deciles("ID.all_variants_vaf.tsv")
  
  colnames(output) <- c("chr", "pos", "ref", "alt", "vaf", "decile")

  # Save the output to a new file in the same directory
  output_file <- "ID.deciles.tsv"
  write.table(output, file = output_file, sep = "\t", row.names = FALSE, quote = FALSE)
  
  RCODE
  >>>  

  runtime {
    memory:  "~{jobMemory} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }
  output {
    File deciles = "ID.deciles.txt"
  }
}


task format4PCGR{
  input {
    File vcfin
    String modules = "neopipe/1.0.0 bcftools/1.9"
    Int jobMemory = 6
    Int timeout = 20
  }
  parameter_meta {
  }
  command <<<
  module load neopipe/1.0.0 bcftools/1.9
  
  python3 /.mounts/labs/gsiprojects/gsi/gsiusers/hdriver/Scripting/NeoAntigen/Install_Mugqic/mugqic_tools-2.12.8.tar.gz/python-tools/format2pcgr.py \
        -i ~{vcfin} \
        -o ID.ensemble.somatic.vt.annot.2callers.TEMP.vcf.gz \
        -f 2 \
        -t TUMOR \
        -v somatic \
        > format2pcgr.log 2>&1
  bcftools view -Oz -i 'TDP>=10 && TVAF>=0.05 && NDP>=10 && NVAF<=0.02' \
        -o ID.ensemble.somatic.vt.annot.2callers.vcf.gz \
        ID.ensemble.somatic.vt.annot.2callers.TEMP.vcf.gz
  tabix -pvcf ID.ensemble.somatic.vt.annot.2callers.vcf.gz
  
  >>>
  runtime {
    memory:  "~{jobMemory} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }
  output {
    File vcfout = "ID.ensemble.somatic.vt.annot.2callers.vcf.gz"
    File vcfoutIndex = "ID.ensemble.somatic.vt.annot.2callers.vcf.gz.tbi"
  }
}


task PCGR{
  input {
    File vcf
    File vcfIndex
    String modules = "pcgr/2.0.3"
    Int jobMemory = 6
    Int timeout = 20
  }
  parameter_meta {
  }
  command <<<
  module load pcgr/2.0.3
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
    --sample_id PCGR \
    --vep_dir $VEP_DIR \
    --refdata_dir $REFDATA_DIR \
    --pcgrr_conda $PCGR_ROOT \
    --no_reporting
  
  R --vanilla<<RCODE
    library(pcgrr)
    yaml_fname<-"pcgr/PCGR.pcgr.grch38.conf.yaml"
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
  python3 <<CODE
  ### filter the snv_indel tsv file
  import csv
  import gzip
  import re
  fout = open("pcgr_filter_sites.txt","w")
  with gzip.open("pcgr/PCGR.pcgr.grch38.snv_indel_ann.tsv.gz",'rt') as f:
    reader=csv.DictReader(f,delimiter="\t")
    for row in reader:
      if row['EXONIC_STATUS']=='exonic' or (row['EXONIC_STATUS']=='non-exonic' and row['ACTIONABILITY_TIER']<=2):
        m=re.search('(.+):g\.(\d+)',row['GENOMIC_CHANGE'])
        chrom="chr" + m[1]
        pos=m[2]
        fout.write("%s\t%s\n" %(chrom,pos))
  fout.close()
  CODE
  
  ### now filter based on sites
  bcftools view -R pcgr_filter_sites.txt -o ID.candidate_sites.vcf ~{vcf}


  	
  >>>
  runtime {
    memory:  "~{jobMemory} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }
  output {
    File candidateCalls = "ID.candidate_sites.vcf"
  }
}










