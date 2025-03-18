import sys
import pandas as pd

def extract_hlas(hlafiles, hlacallers):
    hlas=[]
    for file in hlafiles:
        caller=hlacallers.pop(0)
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
    with open("hlastring.txt","w") as hlaout:
        hlaout.write(hlastring)
    hlaout.close()

if __name__ == "__main__":
    hlafiles = sys.argv[1].split()
    hlacallers = sys.argv[2].split()
    extract_hlas(hlafiles, hlacallers)