#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=20000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o slurm_out/MOLM13_CTCF.out 

union_file="/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f2_union_CTCFs/fz_union_combination/f2_union_binding_data/union_binding_occupancy_score_GT3.bed"

prename="MOLM13_CTCF"
bam_file="/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/22_SWI_SNF/f0_data_process/ChIP_seq_GSE113040_process/process_qc_out_chip/MOLM13_CTCF_DMSO/MOLM13_CTCF_DMSO_MOLM13_CTCF_DMSO_ChIP-Seq.bam"

python get_pattern_near_site_readcount_write_out.py -s hg38 -f bam -m -w 1000 -b 200 -i ${union_file} -t ${bam_file} -o f1_chip_binding_pattern/${prename}_binding_pattern_w1000.csv
python get_pattern_near_site_readcount_write_out.py -s hg38 -f bam -m -w 2000 -b 200 -i ${union_file} -t ${bam_file} -o f1_chip_binding_pattern/${prename}_binding_pattern_w2000.csv
