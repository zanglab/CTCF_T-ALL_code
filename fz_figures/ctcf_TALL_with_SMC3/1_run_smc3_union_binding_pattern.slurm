#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=60000
#SBATCH -t 12:00:00
#SBATCH -p gpu
#SBATCH -A zanglab
#SBATCH -o out1.log

union_binding="/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f2_union_CTCFs/fz_union_combination/f2_union_binding_data/union_binding_occupancy_score_GT3_mid_position.bed"
smc3_results_dir="/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/9_get_CTCF_signals_plus_appended_data/f9_meta_analysis/f3_TF_co-binding/f2_SMC3/f1_combine_2_replicates/chilin_out_combined_replicates/"

mkdir f1_binding_pattern_csv

python get_pattern_near_site_readcount_write_out.py -s hg38 -f bam -m -w 1000 -b 200 -i ${union_binding} -t ${smc3_results_dir}/attic/SMC3_combined_treatment.bam -o f1_binding_pattern_csv/SMC3_CUTLL1_on_union.csv


