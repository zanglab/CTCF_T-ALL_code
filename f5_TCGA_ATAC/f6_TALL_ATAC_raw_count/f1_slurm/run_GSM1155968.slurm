#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=20000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o f1_slurm/slurm_GSM1155968.out

python get_raw_count_on_regions_readcount.py -s hg38 -f bed -g 0 -i /nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f2_union_CTCFs/fz_union_combination/f2_union_binding_data/union_binding_occupancy_score_GT3.bed -t /nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/12_new_gained_lost_from_all_features/f13_ATAC_Jurkat_CD4/f1_proprocess/GSM1155968/GSM1155968_CD4+_ATACseq_Day3_Rep1_PEnoM.bed -o f1_TALL_ATAC_raw_count_on_union/GSM1155968.txt