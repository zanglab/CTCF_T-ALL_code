#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=20000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab


cat f1_peak2000_datasets_union_summits/collection_all_CTCF_narrowPeak.bed |awk '{printf"%s\t%d\t%d\t%s\t%.2f\t%s\t%.2f\n",$1,$2+$10,$2+$10+1,"N",$7,".",$9}' > f1_peak2000_datasets_union_summits/collection_all_CTCF_summits.bed
