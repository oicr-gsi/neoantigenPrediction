import pandas as pd
import re
import sys

# sys.argv[1] = space-separated hlafiles
hlafiles = sys.argv[1].split()
# sys.argv[2] = space-separated hlacallers
hlacallers = sys.argv[2].split()
# sys.argv[3] = output file name (with .txt suffix)
outputFile = sys.argv[3]

# Regex patterns for allele validation
# Matches HLA-A*02:01 or HLA-B*07:02:01
t1k_pattern = re.compile(r'^HLA-[ABC]\*\d{2}:\d{2}(:\d{2,3}){0,2}[A-Z]*$')
# Matches A*02:01 or B*07:02:01
optitype_pattern = re.compile(r'^[ABC]\*\d{2}:\d{2}(:\d{2,3}){0,2}[A-Z]*$')

# this code will extract top HLAs from the various outputs, and store in a common format
all_hlas = []

for file, caller in zip(hlafiles, hlacallers):

    if caller == "t1k":
        # For t1k output, use only the top 3 lines (corresponding to HLA-A, -B, -C)
        try:
            with open(file) as f:
                lines = f.readlines()[0:3]
                for line in lines:
                    fields = line.strip().split()
                    if len(fields) < 6:
                        continue

                    # Only process genes A, B, or C
                    gene_name = fields[0]
                    if gene_name not in ["HLA-A", "HLA-B", "HLA-C"]:
                        continue

                    # Extract both alleles from the line and validate
                    for idx in [2, 5]:
                        allele_raw = fields[idx]
                        if t1k_pattern.match(allele_raw):
                            parts = allele_raw.split("*")
                            gene = parts[0].replace("HLA-", "")
                            allele = ":".join(parts[1].split(":")[0:2])
                            all_hlas.append([gene, f"HLA-{gene}{allele}"])
                        else:
                            print(f"Invalid allele format: {allele_raw}")
        except Exception as e:
            print(f"Error reading T1K file {file}: {e}")

    elif caller == "optitype":
        # For OptiType output, parse the second line which contains the alleles
        try:
            with open(file) as f:
                lines = f.readlines()
                if len(lines) < 2:
                    continue
                fields = lines[1].strip().split()

                # Define index positions for A, B, and C alleles in the optitype output
                gene_cols = [('A', 1, 2), ('B', 3, 4), ('C', 5, 6)]
                for gene, idx1, idx2 in gene_cols:
                    if len(fields) <= max(idx1, idx2):
                        continue
                    
                    # Validate and extract allele1
                    for allele_raw in [fields[idx1], fields[idx2]]:
                        if optitype_pattern.match(allele_raw):
                            allele = allele_raw.split("*")[1]
                            all_hlas.append([gene, f"HLA-{gene}{allele}"])
                        else:
                            print(f"Invalid optitype allele: {allele_raw}")
        except Exception as e:
            print(f"Error reading OptiType file {file}: {e}")

    else:
        print(f"Unknown caller: {caller}")
        continue

# Convert list of alleles into a DataFrame
df = pd.DataFrame(all_hlas, columns=['Gene', 'HLA'])

# Count the number of times each allele appears per gene
dfcounts = df.groupby(['Gene', 'HLA']).size().reset_index(name='Count')

# Sort by gene, then by count (descending), then by allele name (ascending)
dfcounts = dfcounts.sort_values(['Gene', 'Count', 'HLA'], ascending=[True, False, True])

# Select top 2 HLA alleles for each gene (A, B, C)
final_hlas = []

# Get top 2 alleles for this gene, if only one allele is found, duplicate it
for gene in dfcounts["Gene"].unique():
    top_alleles = dfcounts[dfcounts["Gene"] == gene][:2]["HLA"].tolist()
    if len(top_alleles) == 1:
        top_alleles.append(top_alleles[0])
    final_hlas.extend(top_alleles)

# Create final HLA string and write to file
hlastring = " ".join(final_hlas)
try:
    with open(outputFile, "w") as hlaout:
        hlaout.write(hlastring)
    print(f"File {outputFile} has been written")
except Exception as e:
    print(f"Error writing output file: {e}")