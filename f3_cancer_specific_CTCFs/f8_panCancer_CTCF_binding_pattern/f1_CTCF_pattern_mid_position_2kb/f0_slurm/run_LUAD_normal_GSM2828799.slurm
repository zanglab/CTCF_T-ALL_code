#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=20000
#SBATCH -t 12:00:00
#SBATCH -p gpu
#SBATCH -A zanglab
#SBATCH -o f0_slurm_append/slurm_LUAD_normal_GSM2828799.out

python get_pattern_near_site_readcount_write_out.py -s hg38 -f bam -m -w 1000 -b 200 -i union_binding_occupancy_score_GT3_mid_position.bed -t /nv/vol190/zanglab/zw5j/projects_data/T_ALL_CTCF_panos/CTCF/bam/GSM2828799.bam -o f1_ctcf_binding_pattern_csv/LUAD_normal_GSM2828799.txt