#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=40000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o f1_slurm/slurm_JURKAT_H3K27ac_53708.out

python get_RPKM_on_regions_readcount.py -s hg38 -f bam -m -e 10000 -i /nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f2_union_CTCFs/fz_union_combination/f2_union_binding_data/union_binding_occupancy_score_GT3_mid_position.bed -t /nv/vol190/zanglab/cz3d/cistrome_data/bam/53708_treat_rep1.bam -o /nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f10_TF_binding/f3_TALL_Notch/f2_notch_H3K27ac/f1_HM_RPKM/JURKAT_H3K27ac_53708.txt