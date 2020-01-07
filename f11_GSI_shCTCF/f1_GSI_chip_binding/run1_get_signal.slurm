#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=40000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab

#Run program

# union_bed="/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/9_get_CTCF_signals_plus_appended_data/f4_lost_on_consistent_binding/f1_consistent_bindings/f0_infiles/union_summits_width_150.bed"
union_bed="/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f2_union_CTCFs/fz_union_combination/f2_union_binding_data/union_binding_occupancy_score_GT3.bed"

dmso_bam="/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/12_new_gained_lost_from_all_features/f6_Jurkat_GSI/f0_QC/Jurkat_GSI_chilin_out/JURKAT_DMSO/attic/JURKAT_DMSO_treatment.bam"
gsi_3d_bam="/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/12_new_gained_lost_from_all_features/f6_Jurkat_GSI/f0_QC/Jurkat_GSI_chilin_out/JURKAT_GSI_3d/attic/JURKAT_GSI_3d_treatment.bam"
gsi_3d_w4h_bam="/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/12_new_gained_lost_from_all_features/f6_Jurkat_GSI/f0_QC/Jurkat_GSI_chilin_out/JURKAT_GSI_3d_w4hr/attic/JURKAT_GSI_3d_w4hr_treatment.bam"

mkdir f1_signal_csv

time python get_RPKM_on_regions_readcount.py -i $union_bed -t $dmso_bam -s hg38 -f bam -o f1_signal_csv/jurkat_dmso.csv
time python get_RPKM_on_regions_readcount.py -i $union_bed -t $gsi_3d_bam -s hg38 -f bam -o f1_signal_csv/jurkat_gsi_3d.csv
time python get_RPKM_on_regions_readcount.py -i $union_bed -t $gsi_3d_w4h_bam -s hg38 -f bam -o f1_signal_csv/jurkat_gsi_3d_w4h.csv

