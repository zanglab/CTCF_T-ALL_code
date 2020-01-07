#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=40000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o py1.out
#Run program

# union_bed="/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/9_get_CTCF_signals_plus_appended_data/f4_lost_on_consistent_binding/f1_consistent_bindings/f0_infiles/union_summits_width_150.bed"
union_bed="/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f2_union_CTCFs/fz_union_combination/f2_union_binding_data/union_binding_occupancy_score_GT3.bed"
patient_summit_dir="/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/12_new_gained_lost_from_all_features/f4_patient_peak_overlap/patient_summits/"
time python find_overlap_keep_info_NOT_sep_strand_lastColMarked.py -a $union_bed -b ${patient_summit_dir}/PD31_summits.bed -s hg38 -p patient_overlap_out/PD31_overlap.bed -e2 75
# time python find_overlap_keep_info_NOT_sep_strand_lastColMarked.py -a $union_bed -b ${patient_summit_dir}/PD40_summits.bed -s hg38 -p patient_overlap_out/PD40_overlap.bed -e2 75
time python find_overlap_keep_info_NOT_sep_strand_lastColMarked.py -a $union_bed -b ${patient_summit_dir}/PD9_summits.bed  -s hg38 -p patient_overlap_out/PD9_overlap.bed -e2 75
# time python find_overlap_keep_info_NOT_sep_strand_lastColMarked.py -a $union_bed -b ${patient_summit_dir}/PTBG_summits.bed -s hg38 -p patient_overlap_out/PTBG_overlap.bed -e2 75
