import sys
import csv
import gzip
import re

# outputFilePrefix
outputFilePrefix = sys.argv[1]

# Open output file to write filtered sites
with open(f"{outputFilePrefix}.pcgr_filter_sites.txt", "w") as fout:
    with gzip.open(f"pcgr/{outputFilePrefix}.pcgr.grch38.snv_indel_ann.tsv.gz", 'rt') as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            if row['EXONIC_STATUS'] == 'exonic' or (row['EXONIC_STATUS'] == 'non-exonic' and int(row['ACTIONABILITY_TIER']) <= 2):
                # Extract chromosome and position from the GENOMIC_CHANGE field
                m = re.search('(.+):g\.(\d+)', row['GENOMIC_CHANGE'])
                chrom = "chr" + m[1]
                pos = m[2]
                fout.write(f"{chrom}\t{pos}\n")
