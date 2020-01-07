#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=100000
#SBATCH -t 24:00:00
#SBATCH -p largemem
#SBATCH -A zanglab
#SBATCH -o get_pattern.out 

mkdir f1_chip_binding_pattern

union_file="/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f2_union_CTCFs/fz_union_combination/f2_union_binding_data/union_binding_occupancy_score_GT3.bed"

dmso_bam="/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/12_new_gained_lost_from_all_features/f6_Jurkat_GSI/f0_QC/Jurkat_GSI_chilin_out/JURKAT_DMSO/attic/JURKAT_DMSO_treatment.bam"
gsi_3d_bam="/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/12_new_gained_lost_from_all_features/f6_Jurkat_GSI/f0_QC/Jurkat_GSI_chilin_out/JURKAT_GSI_3d/attic/JURKAT_GSI_3d_treatment.bam"
gsi_3d_w4h_bam="/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/12_new_gained_lost_from_all_features/f6_Jurkat_GSI/f0_QC/Jurkat_GSI_chilin_out/JURKAT_GSI_3d_w4hr/attic/JURKAT_GSI_3d_w4hr_treatment.bam"


prename="dmso"
python get_pattern_near_site_readcount_write_out.py -s hg38 -f bam -m -w 1000 -b 200 -i ${union_file} -t ${dmso_bam} -o f1_chip_binding_pattern/${prename}_binding_pattern_w1000.csv

prename="gsi_3d"
python get_pattern_near_site_readcount_write_out.py -s hg38 -f bam -m -w 1000 -b 200 -i ${union_file} -t ${gsi_3d_bam} -o f1_chip_binding_pattern/${prename}_binding_pattern_w1000.csv

prename="gsi_3d_w4h"
python get_pattern_near_site_readcount_write_out.py -s hg38 -f bam -m -w 1000 -b 200 -i ${union_file} -t ${gsi_3d_w4h_bam} -o f1_chip_binding_pattern/${prename}_binding_pattern_w1000.csv

