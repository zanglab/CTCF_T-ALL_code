#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=100000
#SBATCH -t 24:00:00
#SBATCH -p largemem
#SBATCH -A zanglab 
#SBATCH -o out_run5.log


ctcf_bam_file="/nv/vol190/zanglab/zw5j/projects_data/T_ALL_CTCF_panos/CTCF/TALL_bam_merge/Jurkat_merge.bam"
atac_bed_file="data/Jurkat_ATAC_PEnoM.bed"
# region_file_mid="union_CTCF_No_Jurkat_DNAmethylation_with_motif_ATAC_peak_overlapped_mid_position.bed"
# region_file="union_CTCF_No_Jurkat_DNAmethylation_with_motif_ATAC_peak_overlapped.bed"
union_binding_mid="/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f2_union_CTCFs/fz_union_combination/f2_union_binding_data/union_binding_occupancy_score_GT3_mid_position.bed"
union_binding="/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f2_union_CTCFs/fz_union_combination/f2_union_binding_data/union_binding_occupancy_score_GT3.bed"

# python get_RPKM_on_regions_readcount.py -i $union_binding_mid -t $ctcf_bam_file -s hg38 -f bam -m -e 200 -o binding_csv/ctcf_RPKM_e200.csv &
# python get_RPKM_on_regions_readcount.py -i $union_binding_mid -t $atac_bed_file -s hg38 -f bed -m -e 200 -o binding_csv/atac_RPKM_e200.csv
# wait
# python get_RPKM_on_regions_readcount.py -i $union_binding -t $ctcf_bam_file -s hg38 -f bam -o binding_csv/ctcf_RPKM_union_site.csv &
# python get_RPKM_on_regions_readcount.py -i $union_binding -t $atac_bed_file -s hg38 -f bed -o binding_csv/atac_RPKM_union_site.csv

python get_pattern_near_site_readcount_write_out.py -s hg38 -f bam -m -w 1000 -b 200 -i ${union_binding_mid} -t ${ctcf_bam_file} -o binding_csv/ctcf_binding_pattern_w1000.csv
python get_pattern_near_site_readcount_write_out.py -s hg38 -f bed -m -w 1000 -b 200 -i ${union_binding_mid} -t ${atac_bed_file} -o binding_csv/atac_binding_pattern_w1000.csv
