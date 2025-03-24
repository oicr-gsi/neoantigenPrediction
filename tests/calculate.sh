#!/bin/bash
set -o nounset
set -o errexit
set -o pipefail

cd $1

### exclude columns from tsv with precision values
cat *.tsv | cut -f 1-10,13-26,28,30-36 | md5sum

### below will check the xlsx, but will not run due to precision values existing in worksheets which is hard to manage
#xlsxfile=`ls *.xlsx | head -n 1`
#unzip -q $xlsxfile -d xlsx
#find xlsx -name *.xml | xargs md5sum | grep -v core.xml

