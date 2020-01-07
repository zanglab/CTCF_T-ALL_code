#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=40000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab

#Run program

union_binding="/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f2_union_CTCFs/fz_union_combination/f2_union_binding_data/union_binding_occupancy_score_GT3.bed"
union_binding_mid_pos="/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f2_union_CTCFs/fz_union_combination/f2_union_binding_data/union_binding_occupancy_score_GT3_mid_position.bed"

aml_smarca4_bam="/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/22_SWI_SNF/f0_data_process/ChIP_seq_GSE113040_process/process_qc_out_chip/EOL1_SMARCA4_Naive/EOL1_SMARCA4_Naive_EOL1_SMARCA4_Naive_ChIP-Seq.bam"
time python get_RPKM_on_regions_readcount.py -i $union_binding -t $aml_smarca4_bam -o binding_RPKM/EOL1_SMARCA4_binding_region_RPKM.csv -s hg38 -f bam 
time python get_RPKM_on_regions_readcount.py -i $union_binding_mid_pos -t $aml_smarca4_bam -o binding_RPKM/EOL1_SMARCA4_mid_500bp_es_RPKM.csv -s hg38 -f bam -m -e 500

hct116_smarca4_bam="/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/23_HCT116_SMARCA4/f0_data_process/process_qc_out/HCT116_SMARCA4_parental/HCT116_SMARCA4_parental_HCT116_SMARCA4_parental_GSM1835989.bam"
time python get_RPKM_on_regions_readcount.py -i $union_binding -t $hct116_smarca4_bam -o binding_RPKM/HCT116_SMARCA4_binding_region_RPKM.csv -s hg38 -f bam 
time python get_RPKM_on_regions_readcount.py -i $union_binding_mid_pos -t $hct116_smarca4_bam -o binding_RPKM/HCT116_SMARCA4_mid_500bp_es_RPKM.csv -s hg38 -f bam -m -e 500

