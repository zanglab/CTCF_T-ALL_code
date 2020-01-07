#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=60000
#SBATCH -t 24:00:00
#SBATCH -p largemem
#SBATCH -A zanglab 
#SBATCH -o out_run4.log


union_binding_mid="/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f2_union_CTCFs/fz_union_combination/f2_union_binding_data/union_binding_occupancy_score_GT3_mid_position.bed"
union_binding="/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f2_union_CTCFs/fz_union_combination/f2_union_binding_data/union_binding_occupancy_score_GT3.bed"

atac_peak="/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f11_GSI_shCTCF/f5_GSI_ATAC_chip_20191029/f0_process/GSI_ATAC_peak.narrowPeak"

python find_overlap_keep_info_NOT_sep_strand.py -a ${union_binding_mid} -b ${atac_peak} -s hg38 -p union_binding_mid_GSI_ATAC_peak_overlapped.bed -e2 1000
python find_overlap_keep_info_NOT_sep_strand.py -a ${union_binding} -b ${atac_peak} -s hg38 -p union_binding_GSI_ATAC_peak_overlapped.bed -e2 1000

tall_gained="/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f3_cancer_specific_CTCFs/f3_cancer_specific_gained_lost/f1_cancer_specific_binding/T-ALL_gained.bed"
python find_overlap_keep_info_NOT_sep_strand.py -a ${tall_gained} -b union_binding_mid_GSI_ATAC_peak_overlapped.bed -s hg38 -p gained_overlapped.bed