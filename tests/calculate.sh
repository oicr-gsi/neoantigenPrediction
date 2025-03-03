#!/bin/bash
set -o nounset
set -o errexit
set -o pipefail

cd $1

md5sum MoHQ-CM-1-180_neoantigenPredictions.tsv
md5sum MoHQ-CM-1-180_neoantigenPredictions.xlsx

