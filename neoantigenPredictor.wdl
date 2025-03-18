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
      hlafiles = HLAFiles["files"],
      hlacallers = HLAFiles["callers"]
  }
  
  ### prepare for PCGR by adding in required INFO fields : TDP,TVAF,NDP,NVAF
  call format2pcgr {
    input:
      vcfin = DNAVariantCalls["vcf"],
      tumorId = DNAVariantCalls["tumorId"]
  }


  ### generate Deciles from the formatted vcf file
  call vafDeciles {
    input: 
      vcf= format2pcgr.vcfout
  }

  ### call the PCGR software to determine which variants to keep as candidate sites
  call PCGR {
    input:
      vcf = format2pcgr.vcfout,
      vcfIndex = format2pcgr.vcfoutIndex
 }


 call vepAnnotate {
    input:
      vcf = PCGR.candidateCalls,
      refFasta = resources[reference].refFasta
  }

  call getPeptides {
    input: 
      vcf = vepAnnotate.annotatedCandidateCalls,
      tumorId = DNAVariantCalls["tumorId"]

  }

  call formatCalls {
    input:
      peptides = getPeptides.peptides,
      vcf = vepAnnotate.annotatedCandidateCalls
  }

  call ExpressionDeciles {
    input:
      tsv = RNAAbundance
  }
  
  call rnaseqVariants {
    input:
      vcf = RNAVariantCalls["vcf"]
  }
  
  call mergePredictorInputs {
     input:
       variants_peptides = formatCalls.tsv,
       variant_deciles = vafDeciles.deciles,
       expression_deciles = ExpressionDeciles.deciles,
       rnaseq_variants = rnaseqVariants.tsv
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



task mergePredictorOutputs {
   input {
     Array[File] predictorOutputs
     File predictorInputTSV
     String outputFilePrefix
     String modules = "neopipe/1.0.0 sb-neoantigen-models/1.0.0"
     Int jobMemory = 6
     Int timeout = 20	
   }

   command<<<
   python scripts/merge_predictor_outputs.py "~{sep=" " predictorOutputs}" "~{predictorInputTSV}" "~{outputFilePrefix}"
   >>>
}  

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


task chunkPredictorInputFile {
   input{
     File xls
     Int chunksize
     String modules = "neopipe/1.0.0 sb-neoantigen-models/1.0.0"
     Int jobMemory = 6
     Int timeout = 20	
   }

  command<<<
  python3 scripts/chunkPredictorInputFile.py "~{xls}" "~{chunksize}"
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


task extractHLAs{
  input{
    Array[File] hlafiles
    Array[String] hlacallers
    String modules = "neopipe/1.0.0" 
    Int jobMemory = 6
    Int timeout = 20	
  }
  command<<<
  python3 scripts/extractHLAs.py "~{sep=" " hlafiles}" "~{sep=" " hlacallers}"
  >>>

  runtime {
    memory:  "~{jobMemory} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }
  
  output {
    String hlas = read_string("hlastring.txt")
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


task mergePredictorInputs{
   input {
     File variants_peptides
     File variant_deciles
     File expression_deciles
     File rnaseq_variants
     String modules = "neopipe/1.0.0"
     Int jobMemory = 6
     Int timeout = 20	  
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
   write.table(m, "ID.output.merged.tsv", quote=F, row.names=F, sep="\t")
   # Format for .xlsx input
   df=data.frame(paste(m\$chr,m\$pos,m\$ref,m\$alt,sep=";"),m\$wt_peptides,m\$mt_peptides,m\$vaf_decile,m\$expression_decile,m\$variant_in_rna)
   colnames(df)=c("Unique identifier","Wt nmer","Mut nmer","Exome VAF decile","Gene expression decile","Present in RNA-seq data")
   write_xlsx(df, "ID.output.merged.xls")
   RCODE
   >>>  

  runtime {
    memory:  "~{jobMemory} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }
  output {
    File tsv = "ID.output.merged.tsv"
    File xls = "ID.output.merged.xls"
  }

}


task rnaseqVariants{
  input {
    File vcf
    String modules = "bcftools/1.9"
    Int jobMemory = 6
    Int timeout = 20
  }
  command<<<
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t1\n' ~{vcf} > ID.rnaseq_variants.tsv
    >>>  

  runtime {
    memory:  "~{jobMemory} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }
  output {
    File tsv = "ID.rnaseq_variants.tsv"
  }
}


task formatCalls{
  input {
    File peptides
    File vcf
    String modules = "bcftools/1.9"
    Int jobMemory = 6
    Int timeout = 20
  }
  

  command<<<
    python3 <<CODE
    import re
    import os

    fout = open("ID.candidateCalls.tsv","w")
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
    File tsv = "ID.candidateCalls.tsv"
  }
}



task ExpressionDeciles{
  input {
    File tsv
    String modules = "neopipe/1.0.0"
    Int jobMemory = 6
    Int timeout = 20
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
  tx2gene_file<- "/.mounts/labs/gsi/src/neoAntigen/Homo_sapiens.GRCh38.Ensembl104.tx2gene"
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
  output_fname<-"ID.expression_deciles_kallisto50_ensembl104.tsv"
  write.table(output, file = output_fname, sep = "\t", row.names = FALSE, quote=F)
  
  RCODE
  >>>  

  runtime {
    memory:  "~{jobMemory} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }
  output {
    File deciles = "ID.expression_deciles_kallisto50_ensembl104.tsv"
  }
}



task getPeptides {
  input{
    File vcf
    String tumorId
    String modules = "pvactools/4.3.0"
    Int jobMemory = 6
    Int timeout = 20
   }
   
   command<<<

   pvacseq generate_protein_fasta \
   -s ~{tumorId} \
   -d 12 \
   ~{vcf} 12 ID.peptides.fa
   
   ## convert to single line per site with WTid,WTseq,MTid,MTseq
   cat ID.peptides.fa | paste - - - - > ID.peptides.0.txt
   ## keep only records where MT != WT peptide
   cat ID.peptides.0.txt | awk '{if($2 != $4) print}' > ID.peptides.1.txt
   ## order by numeric identifier
   cat ID.peptides.1.txt | sort -t . -k2 -g > ID.peptides.2.txt
   ## keep only if unique and peptide length == 25
   cat ID.peptides.2.txt | awk -F'\t' '!seen[$2,$4]++ && length($2) >= 25 && length($4) >= 25' > ID.peptides.txt


   >>>

  runtime {
    memory:  "~{jobMemory} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }
  output {
    File peptides = "ID.peptides.txt"
  }
}

task vepAnnotate{
  input{
    File vcf
    String refFasta
    String modules = "vep/112.0 pcgr/2.0.3 hg38/p12"
    Int jobMemory = 6
    Int timeout = 20
   }
   
   
   String plugins="/.mounts/labs/gsi/src/pvactools/plugins/"
   command<<<
   
   vep \
   --input_file ~{vcf} \
   --output_file ID.vep.vcf \
   --format vcf --vcf --symbol --terms SO --tsl \
   --hgvs --fasta ~{refFasta} \
   --offline \
   --cache --dir_cache $VEP_DIR --cache_version 112 \
   --plugin Frameshift --plugin Wildtype \
   --dir_plugins ~{plugins} --pick --transcript_version \
   --force_overwrite
   >>>

  runtime {
    memory:  "~{jobMemory} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }
  output {
    File annotatedCandidateCalls = "ID.vep.vcf"
  }
}

task vafDeciles{
  input {
    File vcf
    String modules = "neopipe/1.0.0 bcftools/1.9"
    Int jobMemory = 6
    Int timeout = 20
  }
  
  command<<<

  ## extract relevant fields from the vcf records
  bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%TVAF\n' ~{vcf} > ID.all_variants_vaf.tsv
  
  ## generated deciles with R code
  R --vanilla <<RCODE
  library(dplyr)
  
  fin="ID.all_variants_vaf.tsv"
  data <- read.delim(fin, header = FALSE, sep = "\t", col.names = c("chr", "pos", "ref", "alt", "vaf"))

  # Standardize the vaf column
  standardized_col <- scale(data\$vaf)
  
  # Calculate deciles
  deciles <- quantile(standardized_col, probs = seq(0, 1, 0.1))

  # Output results
  output <- data.frame(data[, 1:5], decile = cut(standardized_col, breaks = deciles, labels = FALSE, include.lowest = TRUE))
  
  # Save the output to a new file in the same directory
  fout <- "ID.deciles.tsv"
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
    File deciles = "ID.deciles.tsv"
  }
}


task format2pcgr{
  input {
    File vcfin
    String tumorId
    String modules = "neopipe/1.0.0 bcftools/1.9"
    Int jobMemory = 6
    Int timeout = 20
  }
  parameter_meta {
  }
  command <<<
  
  python3 /.mounts/labs/gsiprojects/gsi/gsiusers/hdriver/Scripting/NeoAntigen/Install_Mugqic/mugqic_tools-2.12.8.tar.gz/python-tools/format2pcgr.py \
        -i ~{vcfin} \
        -o ID.ensemble.somatic.vt.annot.2callers.TEMP.vcf.gz \
        -f 2 \
        -t ~{tumorId} \
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










