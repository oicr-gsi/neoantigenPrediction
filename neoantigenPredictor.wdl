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
    String modules = "neopipe/1.0.0" 
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
  python3 <<CODE
  import pandas as pd
  
  files="~{sep=" " hlafiles}".split()
  callers="~{sep=" " hlacallers}".split()
  
  ### this code will extract top HLAs from the various outputs, and store in a common format
  hlas=[]
  for file in files:
    caller=callers.pop(0)
    if caller == "t1k":
      ### use the top 3 alleles (A,B,C), and show the HLAGene Family (col1, 5th character, A,B or C), and the first and second HLA Allele names (col 3,6)
      with open(file) as f:
        lines=f.readlines()[0:3]
        for line in lines:
          fields=line.replace("*","").split()
          ## get the gene HLA-GENE* as a single character, A,B or C, a class I HLA
          hla_gene=fields[0][4:]
    
          ## limit HLA codes to Field 1 and 2, removing everything after (split on :, then join the first two fields)
          allele1=":".join(fields[2].split(":")[0:2])
          allele2=":".join(fields[5].split(":")[0:2])
          hlas.append([hla_gene,allele1])
          hlas.append([hla_gene,allele2])

    elif caller == "optitype":
      ### A,B and C Alleles are all on one line in columns 2-7 (A1,A2,B1,B2,C1,C2)
      with open(file) as f:
        lines=f.readlines()
        fields=lines[1].replace("*","").split()
        hlas.append(["A","HLA-" + fields[1]])
        hlas.append(["A","HLA-" + fields[2]])
        hlas.append(["B","HLA-" + fields[3]])
        hlas.append(["B","HLA-" + fields[4]])
        hlas.append(["C","HLA-" + fields[5]])
        hlas.append(["C","HLA-" + fields[6]])
    else:
      print("unknown caller " + caller)
      quit()

  ### convert to data frame and get counts and sort by Count, then alphabetically by the HLA
  df = pd.DataFrame(hlas,columns=['Gene','HLA'])
  dfcounts = df.groupby(['Gene','HLA']).size().reset_index(name='Count')
  dfcounts = dfcounts.sort_values(['Gene', 'Count', "HLA"], ascending = [True, False,True])
  
  ### collect final selection to list
  hlas=[]
  for gene in dfcounts["Gene"].unique():
    ## get the top 2 Genes, if only one has been identified, use it twice in the final list
    dff=dfcounts[dfcounts["Gene"]==gene][0:2]
    
    hla_list=dff["HLA"].values.tolist()
    
    if len(hla_list)<2:
        hla_list.extend(hla_list)
    hlas.extend(hla_list)
  
 
  hlastring= " ".join(hlas)
  with open("~{outputFilePrefix}.hlastring.txt","w") as hlaout:
    hlaout.write(hlastring)
  hlaout.close()
  CODE
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
    String modules = "neopipe/1.0.0 bcftools/1.9"
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
  
  format2pcgr \
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
    String modules = "pcgr/2.0.3"
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
  python3 <<CODE
  ### filter the snv_indel tsv file
  import csv
  import gzip
  import re
  fout = open("~{outputFilePrefix}.pcgr_filter_sites.txt","w")
  with gzip.open("pcgr/~{outputFilePrefix}.pcgr.grch38.snv_indel_ann.tsv.gz",'rt') as f:
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
    String modules = "bcftools/1.9"
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
    python3 <<CODE
    import re
    import os

    fout = open("~{outputFilePrefix}.candidateCalls.tsv","w")
    output_fields=["chr","pos","rs_id","ref","alt","type","callers","tumor_dp","normal_dp","tumor_vaf","normal_vaf","consequence","impact","gene","ensembl_id","transcript_id","biotype","nuc_change","aa_change","cdna_pos","cds_pos","protein_pos","amino_acids","codons","wt_peptides","mt_peptides"]
    print(*output_fields,sep="\t",file=fout)
    
    ## load and catalogue the peptides
    peptides={}
    with open("~{peptides}") as pfile:
      for line in pfile:
        line=line.rstrip("\n")
        WTid,WT_peptide,MTid,MT_peptide=line.split("\t")
        d=dict(zip(["prefix","ordinal","gene","tid0","tid1","type","change"],WTid.split(".")))
        d['wt_peptide']=WT_peptide
        d['mt_peptide']=MT_peptide
        d['transcript_id']=d['tid0'] + "." +d['tid1']
        results=re.search("(\d+)(.*)",d['change'])
        d['amino_acids']=results[2]
        d['protein_pos']=results[1]
        key="|".join([d[x] for x in ['gene','transcript_id','protein_pos','amino_acids']])
        peptides[key]=d

      #### keys/fields of interst
      info_keys=["CHROM","POS","ID","REF","ALT","TYPE","CALLERS","TDP","NDP","TVAF","NVAF"]
      csq_keys=["Consequence","IMPACT","SYMBOL","Gene","Feature","BIOTYPE","HGVSc","HGVSp","cDNA_position","CDS_position","Protein_position","Amino_acids","Codons","wt_peptide","mt_peptide"]

      ### get CSQ fields from the header
      headers=os.popen('bcftools view -h "~{vcf}"').readlines()
      for header in headers:
        if "##INFO" in header:
          if "ID=CSQ" in header:
            CSQ_fields_string=re.search("Format: (.*)\"",header)
            CSQ_fields=CSQ_fields_string[1].split("|")

      ### get required fields from each record
      recs=os.popen('bcftools query -f "%CHROM\t%POS\t%ID\t%REF\t%ALT\t%TYPE\t%CALLERS\t%TDP\t%NDP\t%TVAF\t%NVAF\t%CSQ\n" "~{vcf}"').readlines()
      for rec in recs:
        rec=rec.rstrip("\n")    
        ### get the values for this record and store in info_dict
        v=rec.split("\t")
        info_dict=dict(zip(info_keys + ["CSQ"],v))

        ### get the values for the CSQ part of the record, store in csq_dict
        v=info_dict['CSQ'].split("|")
        csq_dict=dict(zip(CSQ_fields,v))

        ### for a key for the peptide dictionary
        peptide_key="|".join([csq_dict[x] for x in ['SYMBOL','Feature','Protein_position','Amino_acids']])

        ### print out the record only if the key is in the peptide dictionary 
        if peptide_key in peptides:
           csq_dict['wt_peptide']=peptides[peptide_key]['wt_peptide']
           csq_dict['mt_peptide']=peptides[peptide_key]['mt_peptide']
           ### get the required info and csq values (including the peptides), and print to the output file
           info_values=[info_dict[x] for x in info_keys]
           csq_values =[csq_dict[x] for x in csq_keys] 
           all_values=info_values + csq_values
           print(*all_values,sep="\t",file=fout)

      fout.close()
    CODE


    
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
    String modules = "neopipe/1.0.0 bcftools/1.9"
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
  
  ## generated deciles with R code
  R --vanilla <<RCODE
  library(dplyr)
  
  fin="~{outputFilePrefix}.all_variants_vaf.tsv"
  data <- read.delim(fin, header = FALSE, sep = "\t", col.names = c("chr", "pos", "ref", "alt", "vaf"))

  # Standardize the vaf column
  standardized_col <- scale(data\$vaf)
  
  # Calculate deciles
  deciles <- quantile(standardized_col, probs = seq(0, 1, 0.1))

  # Output results
  output <- data.frame(data[, 1:5], decile = cut(standardized_col, breaks = deciles, labels = FALSE, include.lowest = TRUE))
  
  # Save the output to a new file in the same directory
  fout <- "~{outputFilePrefix}.deciles.tsv"
  ### this seens to writes out the data without a headerline...is this necessary, or maybe it is done purposely
  write.table(output, file = fout, sep = "\t", row.names = FALSE, quote = FALSE)
  
  RCODE
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
    String modules = "neopipe/1.0.0 ensembl/104-hg38"
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

  ## convert abundance to gene level
  R --vanilla --args ~{tsv} <<RCODE
  #library(rhdf5)
  library(tximport)
  args <- commandArgs(trailingOnly = TRUE)
  
  ### abundance of transcripts is conferted to gene level, with genes represented by multiple transcripts
  ### see https://bioconductor.org/packages/release/bioc/html/tximport.html
  ### parameters include the expression caller (eg kallisto), which aligns to the expected format
  ### paremeters include a mapping file of transcript IDs to gene IDs (tx2gene)
  tx2gene_file<- "$ENSEMBL_ROOT/Homo_sapiens.GRCh38.Ensembl104.tx2gene"
  tx2gene=read.delim(tx2gene_file, as.is=T)
  df <- as.data.frame(tximport(args[1],type = "kallisto", tx2gene = tx2gene))
  
  #write.table(txi.abundance, "gene_abundance.tsv", sep="\t", col.names=T, row.names=F, quote=F)
  
  ### separate out abundance = 0 genes, do not use for deciles
  df0<-df[df\$abundance ==0,]
  df1<-df[df\$abundance !=0,]
  
  standardized_col <- scale(df1\$abundance)
  deciles <- quantile(standardized_col,probs=seq(0,1,0.1))
  
  ### set quantiles for abundance = 0 genes to 1, apply quantiles to abundance >0 genes
  output0 <- data.frame(gene_id=rownames(df0),abundance_TPM=df0\$abundance,deciles=1)
  output1 <- data.frame(gene_id=rownames(df1),abundance_TPM=df1\$abundance,deciles=cut(standardized_col, breaks = deciles, labels = FALSE, include.lowest = TRUE))
  ## combine output
  output <- rbind(output1,output0)
  output[order(output\$gene_id),]
  output_fname<-"~{outputFilePrefix}.expression_deciles_kallisto50_ensembl104.tsv"
  write.table(output, file = output_fname, sep = "\t", row.names = FALSE, quote=F)
  
  RCODE
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
     String modules = "neopipe/1.0.0"
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

   R --vanilla --args ~{variants_peptides} ~{variant_deciles} ~{expression_deciles} ~{rnaseq_variants} <<RCODE
   library(writexl)
   args <- commandArgs(trailingOnly = TRUE)
   
   vafPeptides=read.csv(args[1], header=T, sep="\t", stringsAsFactors=F)
   vafDeciles=read.csv(args[2], header=T, sep="\t", stringsAsFactors=F)
   colnames(vafDeciles)[ncol(vafDeciles)]="vaf_decile"
   expressionDeciles=read.csv(args[3], header=T, sep="\t", stringsAsFactors=F)
   colnames(expressionDeciles)[ncol(expressionDeciles)]="expression_decile"
   variantsRNASeq=read.csv(args[4], header=F, sep="\t", stringsAsFactors=F)
   colnames(variantsRNASeq)=c("chr","pos","ref","alt","variant_in_rna")

   m=merge(vafPeptides,vafDeciles,by=c("chr","pos","ref","alt"), all.x=TRUE)
   m=merge(m,expressionDeciles,by.x="ensembl_id",by.y="gene_id",all.x=TRUE)
   m$expression_decile[is.na(m$expression_decile)]=1
   m=merge(m,variantsRNASeq, by=c("chr","pos","ref","alt"), all.x=TRUE)
   m$variant_in_rna[is.na(m$variant_in_rna)]=0

   # Write in tsv table
   write.table(m, "~{outputFilePrefix}.output.merged.tsv", quote=F, row.names=F, sep="\t")
   # Format for .xlsx input
   df=data.frame(paste(m\$chr,m\$pos,m\$ref,m\$alt,sep=";"),m\$wt_peptides,m\$mt_peptides,m\$vaf_decile,m\$expression_decile,m\$variant_in_rna)
   colnames(df)=c("Unique identifier","Wt nmer","Mut nmer","Exome VAF decile","Gene expression decile","Present in RNA-seq data")
   write_xlsx(df, "~{outputFilePrefix}.output.merged.xls")
   RCODE
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

  
task mergePredictorOutputs {
   input{
     Array[File] predictorOutputs
     File predictorInputTSV
     String outputFilePrefix
     String modules = "neopipe/1.0.0 sb-neoantigen-models/1.0.0"
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
  python3 <<CODE
  import sys
  import pandas as pd
  import numpy as np

  predictorOutputsString = "~{sep=" " predictorOutputs}"
  predictorOutputs=predictorOutputsString.split()
  data={}
  for ds in ["InputData","NmersScored","MmpsScored"]:
      for predictorOutput in predictorOutputs:
          df=pd.read_excel(predictorOutput,sheet_name=ds)
          if ds not in data:
              data[ds]=df
          else:
              data[ds]=pd.concat([data[ds],df],axis=0)

  writer = pd.ExcelWriter('~{outputFilePrefix}_neoantigenPredictions.xlsx', engine='xlsxwriter')
  data["InputData"].to_excel(writer, sheet_name='InputData', index = False)
  data["NmersScored"].to_excel(writer, sheet_name='NmersScored', index = False)
  data["MmpsScored"].to_excel(writer, sheet_name='MmpsScored', index = False)
  writer.save()
  
  
  #### annotation of the predictorInputTSV with information from the excel file
  ####### NmersScored has 1 record per Unique identifier
  ####### MmpsScores had multiple records per Unique identifer, across a range of MMP model score
  ####### merge the two, matching each of the NmersScored mmp scores to the MMP model score, concatenating values if there are multiple matches
  df_nmers = data["NmersScored"]
  df_mmps = data["MmpsScored"]

  ### to capture the list of lists, for conversion to a dataframe
  list_mmps_match=[]
  ### interate across the nmers row by row
  for index,row in df_nmers.iterrows():
    ### list of values to ad to the list_mmps_match
    match=[]
    uid=row['Unique identifier']    
    match.append(uid)
    
    ## adjust precision
    score1=round(row["mmp score 1"],7)
    score2=round(row["mmp score 2"],7)

    for score in [score1,score2]:
      ## match on uid and score
      mmps_rows = df_mmps[ (df_mmps['Unique identifier'] == uid ) & (round(df_mmps['MMP model score'],7) == score) ]
      ### if multiple rows, concatenate the values
      mutant_peptides=";".join(mmps_rows['Mutant peptide'])
      match.append(mutant_peptides)
      hlas=";".join(mmps_rows['HLA'])
      match.append(hlas)

    ### add to the list
    list_mmps_match.append(match)

  ## convert the list to a dataframe
  df_mmps_match=pd.DataFrame(list_mmps_match, columns=["Unique identifier","mmpScore1Pep","mmpScore1HLA","mmpScore2Pep","mmpScore2HLA"])

  ### merge the new columns to df_nmers
  df_nmers_extended=pd.merge(df_nmers,df_mmps_match,on="Unique identifier")

  ### split the Unique identifier to columns, for the final merge, adjust the dtype for pos for merging
  df_nmers_extended[["chr","pos","ref","alt"]]=df_nmers_extended['Unique identifier'].str.split(";",expand=True,)
  df_nmers_extended=df_nmers_extended.drop(["Unique identifier"],axis=1)
  df_nmers_extended["pos"]=df_nmers_extended["pos"].astype('int64')

  ### load the tsv file
  df_tsv=pd.read_csv("~{predictorInputTSV}",sep="\t")
  ### merge in the nmer information
  df_tsv_extended=pd.merge(df_tsv,df_nmers_extended,on=["chr","pos","ref","alt"])

  df_tsv_extended.to_csv("~{outputFilePrefix}_neoantigenPredictions.tsv",index=False,sep="\t")  
  
  CODE
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
     String modules = "neopipe/1.0.0 sb-neoantigen-models/1.0.0"
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
  
  python3 <<CODE
  import sys
  import pandas as pd
  import numpy as np
  import os


  
  df=pd.read_excel("~{xls}")
  rowcount=len(df)
  #### need to round this appropriate
  chunks=int(np.ceil(rowcount/~{chunksize}))
  dfs=np.array_split(df,chunks)

  chunk=0
  for df in dfs:
      chunk = chunk + 1
      ## this function appears to require xlsx extension, but i need to use xls extension for the predictor
      tmpfn="predict_input" + str(chunk) + ".xlsx"
      df.to_excel(tmpfn,index=False)
      os.rename("predict_input" + str(chunk) + ".xlsx", "predict_input" + str(chunk) + ".xls")
     
  CODE
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
  
  
  











