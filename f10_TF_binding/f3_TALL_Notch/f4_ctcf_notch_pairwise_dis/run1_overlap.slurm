#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=20000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o py1.out
#Run program


# binding_bed="/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/13_interval_distribution/f0_union_binding_files/union_binding_singleton_bi_removed_motif_position.bed"
binding_bed="/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f2_union_CTCFs/fz_union_combination/f2_union_binding_data/union_binding_occupancy_score_GT3_mid_position.bed"

cutll1_myc_peak="/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/18_MYC-ChIP_shCTCF-RNA/MYC_ChIP/f0_processing/callpeak/macs2_out/CUTLL1_MYC_summits.bed"
notch_w4h_peak="/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/12_new_gained_lost_from_all_features/f8_cofactor_T-ALL_Notch/T_ALL/f0_dynamic_notch_co_binding/Cutll1_Notch1_files/CUTLL1_NOTCH1_w4h_sort_summits.bed"
dynamic_notch_peak="/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/12_new_gained_lost_from_all_features/f8_cofactor_T-ALL_Notch/T_ALL/f0_dynamic_notch_co_binding/Cutll1_Notch1_files/GSM1252936_dynamic_NOTCH1_peaks_hg38.bed"

mkdir f1_NOTCH1_MYC_nearest_dis_csv
time python find_overlap_keep_info_NOT_sep_strand_lastColMarked_Distance_keep_revised.py -a $binding_bed -b $cutll1_myc_peak -s hg38 -p f1_NOTCH1_MYC_nearest_dis_csv/union_CTCF_nearest_CUTLL1_MYC.csv
time python find_overlap_keep_info_NOT_sep_strand_lastColMarked_Distance_keep_revised.py -a $binding_bed -b $notch_w4h_peak  -s hg38 -p f1_NOTCH1_MYC_nearest_dis_csv/union_CTCF_nearest_NOTCH1_w4h.csv
time python find_overlap_keep_info_NOT_sep_strand_lastColMarked_Distance_keep_revised.py -a $binding_bed -b $dynamic_notch_peak -s hg38 -p f1_NOTCH1_MYC_nearest_dis_csv/union_CTCF_nearest_NOTCH1_dynamic.csv

