#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=20000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab

cat /nv/vol190/zanglab/zw5j/projects_data/T_ALL_CTCF_panos/CTCF/peaks/*_peaks.narrowPeak > f1_peak2000_datasets_union_summits/collection_all_CTCF_narrowPeak.bed