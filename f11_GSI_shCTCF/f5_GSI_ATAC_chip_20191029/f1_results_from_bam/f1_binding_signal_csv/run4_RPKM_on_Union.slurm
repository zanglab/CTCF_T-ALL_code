#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=60000
#SBATCH -t 24:00:00
#SBATCH -p largemem
#SBATCH -A zanglab 
#SBATCH -o out_run4.log


union_binding_mid="/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f2_union_CTCFs/fz_union_combination/f2_union_binding_data/union_binding_occupancy_score_GT3_mid_position.bed"
union_binding="/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f2_union_CTCFs/fz_union_combination/f2_union_binding_data/union_binding_occupancy_score_GT3.bed"

dmso_reads="/nv/vol190/zanglab/shared/TALL_mayo/Transfer/A1_DMSO_3d_NoDup.bam"
gsi_3d_reads="/nv/vol190/zanglab/shared/TALL_mayo/Transfer/B1_Gammai_3d_NoDup.bam"
gsi_wo_reads="/nv/vol190/zanglab/shared/TALL_mayo/Transfer/C1_Gammai_Wash_NoDup.bam"

mkdir binding_csv

python get_RPKM_on_regions_readcount.py -i $union_binding_mid -t $dmso_reads   -s hg38 -f bam -m -e 150 -o binding_csv/dmso_ATAC_RPKM_e150.csv
python get_RPKM_on_regions_readcount.py -i $union_binding_mid -t $gsi_3d_reads -s hg38 -f bam -m -e 150 -o binding_csv/gsi_3d_ATAC_RPKM_e150.csv
python get_RPKM_on_regions_readcount.py -i $union_binding_mid -t $gsi_wo_reads -s hg38 -f bam -m -e 150 -o binding_csv/gsi_wo_ATAC_RPKM_e150.csv

# python get_RPKM_on_regions_readcount.py -i $union_binding_mid -t $dmso_reads   -s hg38 -f bam -m -e 200 -o binding_csv/dmso_ATAC_RPKM_e200.csv
# python get_RPKM_on_regions_readcount.py -i $union_binding_mid -t $gsi_3d_reads -s hg38 -f bam -m -e 200 -o binding_csv/gsi_3d_ATAC_RPKM_e200.csv
# python get_RPKM_on_regions_readcount.py -i $union_binding_mid -t $gsi_wo_reads -s hg38 -f bam -m -e 200 -o binding_csv/gsi_wo_ATAC_RPKM_e200.csv

# python get_RPKM_on_regions_readcount.py -i $union_binding_mid -t $dmso_reads   -s hg38 -f bam -m -e 500 -o binding_csv/dmso_ATAC_RPKM_e500.csv
# python get_RPKM_on_regions_readcount.py -i $union_binding_mid -t $gsi_3d_reads -s hg38 -f bam -m -e 500 -o binding_csv/gsi_3d_ATAC_RPKM_e500.csv
# python get_RPKM_on_regions_readcount.py -i $union_binding_mid -t $gsi_wo_reads -s hg38 -f bam -m -e 500 -o binding_csv/gsi_wo_ATAC_RPKM_e500.csv
# 
# python get_RPKM_on_regions_readcount.py -i $union_binding -t $dmso_reads   -s hg38 -f bam -o binding_csv/dmso_ATAC_RPKM_binding_site.csv
# python get_RPKM_on_regions_readcount.py -i $union_binding -t $gsi_3d_reads -s hg38 -f bam -o binding_csv/gsi_3d_ATAC_RPKM_binding_site.csv
# python get_RPKM_on_regions_readcount.py -i $union_binding -t $gsi_wo_reads -s hg38 -f bam -o binding_csv/gsi_wo_ATAC_RPKM_binding_site.csv


