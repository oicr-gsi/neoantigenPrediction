#!/bin/bash
set -o nounset
set -o errexit
set -o pipefail

cd $1

md5sum *.tsv

xlsxfile=`ls *.xlsx | head -n 1`
unzip -q $xlsxfile -d xlsx
find xlsx -name *.xml | xargs md5sum | grep -v core.xml

