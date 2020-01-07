#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=100000
#SBATCH -t 24:00:00
#SBATCH -p largemem
#SBATCH -A zanglab
#SBATCH -o slurm_out/TALL_CTCF.out 

# mkdir f1_chip_binding_pattern
# mkdir slurm_out

union_file="/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f2_union_CTCFs/fz_union_combination/f2_union_binding_data/union_binding_occupancy_score_GT3.bed"

prename="cutll1_CTCF"
bam_file="/nv/vol190/zanglab/zw5j/projects_data/T_ALL_CTCF_panos/CTCF/TALL_bam_merge/CUTLL1_merge.bam"
python get_pattern_near_site_readcount_write_out.py -s hg38 -f bam -m -w 1000 -b 200 -i ${union_file} -t ${bam_file} -o f1_chip_binding_pattern/${prename}_binding_pattern_w1000.csv
# python get_pattern_near_site_readcount_write_out.py -s hg38 -f bam -m -w 2000 -b 200 -i ${union_file} -t ${bam_file} -o f1_chip_binding_pattern/${prename}_binding_pattern_w2000.csv


# prename="jurkat_CTCF"
# bam_file="/nv/vol190/zanglab/zw5j/projects_data/T_ALL_CTCF_panos/CTCF/TALL_bam_merge/Jurkat_merge.bam"
# python get_pattern_near_site_readcount_write_out.py -s hg38 -f bam -m -w 1000 -b 200 -i ${union_file} -t ${bam_file} -o f1_chip_binding_pattern/${prename}_binding_pattern_w1000.csv
