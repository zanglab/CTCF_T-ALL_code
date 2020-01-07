#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=40000
#SBATCH -t 24:00:00
#SBATCH -p gpu
#SBATCH -A zanglab 
#SBATCH -o run1_append.out


# ==== PRAD ====

motif_bed_file="/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f2_union_CTCFs/fz_union_combination/f2_union_binding_data/union_binding_occupancy_score_GT3_mid_position.bed"
bdg_dir="/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f8_DNA_methylation/f1_data_process_cov_new/PanCancer"
    

time python get_loci_signal_on_regions_scores.py -i $motif_bed_file -t ${bdg_dir}/PRAD_GEO_GSE86832/GSE86832_bigTable_hg38_LNCaP_cov5_sorted.bdg  -o f1_signal_binding_csv/LNCaP_cov5_DNAme_unionBinding.csv -s hg38 -c 0 -m -e 150
time python get_loci_signal_on_regions_scores.py -i $motif_bed_file -t ${bdg_dir}/prostate_gland_ENCODE/ENCFF027KTR_cov5.bdg  -o f1_signal_binding_csv/prostate_gland_cov5_DNAme_unionBinding.csv -s hg38 -c 0 -m -e 150




