#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=20000
#SBATCH -t 6:00:00
#SBATCH -p standard 
#SBATCH -A zanglab
#SBATCH -o run_py1_slurm/GSM2825675.out

time python count_reads_on_regions.py -i /nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f2_union_CTCFs/f1_union_all_summits/f1_peak2000_datasets_union_summits/union_summits_fe4_width_150_sorted_natural.bed -b /nv/vol190/zanglab/zw5j/projects_data/T_ALL_CTCF_panos/CTCF/bam/GSM2825675.bam -f bam -s hg38 -o f1_RPKM_csv/GSM2825675.csv