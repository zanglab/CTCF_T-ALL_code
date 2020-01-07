#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=40000
#SBATCH -t 24:00:00
#SBATCH -p largemem
#SBATCH -A zanglab
#SBATCH -o f0_slurm/slurm_BRCA_Breast_normal_H3K27me3_breast_epithelium_H3K27me3_ENCFF114URN.out

time python get_pattern_near_site_readcount_write_out.py -s hg38 -f bam -m -w 1000 -b 200 -i /nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f2_union_CTCFs/fz_union_combination/f2_union_binding_data/union_binding_occupancy_score_GT3_mid_position.bed -t /nv/vol190/zanglab/zw5j/projects_data/T_ALL_CTCF_panos/Public/histone_modification/breast_epithelium_H3K27me3_ENCFF114URN.bam -o /nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f9_histone_modification/f1_binding_pattern_on_union/f1_HM_pattern_csv/BRCA_Breast_normal_H3K27me3_breast_epithelium_H3K27me3_ENCFF114URN.txt

