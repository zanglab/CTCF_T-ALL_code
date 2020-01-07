#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=80000
#SBATCH -t 24:00:00
#SBATCH -p largemem
#SBATCH -A zanglab 
#SBATCH -o out_run4.log


union_binding_mid="/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f2_union_CTCFs/fz_union_combination/f2_union_binding_data/union_binding_occupancy_score_GT3_mid_position.bed"
union_binding="/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f2_union_CTCFs/fz_union_combination/f2_union_binding_data/union_binding_occupancy_score_GT3.bed"

union_peak_overlapped="union_binding_GSI_ATAC_peak_overlapped.bed"
union_mid_peak_overlapped="union_binding_mid_GSI_ATAC_peak_overlapped.bed"

mapping_dir="/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f11_GSI_shCTCF/f5_GSI_ATAC_chip_20191029/f0_process/mapping_callpeak"
dmso_reads=${mapping_dir}/"A1_DMSO_3d/A1_DMSO_3d.bam"
gsi_3d_reads=${mapping_dir}/"B1_Gammi_3d/B1_Gammi_3d.bam"
# gsi_3d_reads=${mapping_dir}/"B1_Gammi_3d_saved/B1_Gammi_3d_raw.bam"
gsi_wo_reads=${mapping_dir}/"C1_Gammi_Wash/C1_Gammi_Wash.bam"


mkdir binding_csv

# python get_RPKM_on_regions_readcount.py -i $union_mid_peak_overlapped -t $dmso_reads   -s hg38 -f bam -m -e 200 -o binding_csv/dmso_ATAC_RPKM_peak_overlapped_union_e200.csv
# python get_RPKM_on_regions_readcount.py -i $union_mid_peak_overlapped -t $gsi_3d_reads -s hg38 -f bam -m -e 200 -o binding_csv/gsi_3d_ATAC_RPKM_peak_overlapped_union_e200.csv
# python get_RPKM_on_regions_readcount.py -i $union_mid_peak_overlapped -t $gsi_wo_reads -s hg38 -f bam -m -e 200 -o binding_csv/gsi_wo_ATAC_RPKM_peak_overlapped_union_e200.csv
# 
# python get_RPKM_on_regions_readcount.py -i $union_mid_peak_overlapped -t $dmso_reads   -s hg38 -f bam -m -e 500 -o binding_csv/dmso_ATAC_RPKM_peak_overlapped_union_e500.csv
# python get_RPKM_on_regions_readcount.py -i $union_mid_peak_overlapped -t $gsi_3d_reads -s hg38 -f bam -m -e 500 -o binding_csv/gsi_3d_ATAC_RPKM_peak_overlapped_union_e500.csv
# python get_RPKM_on_regions_readcount.py -i $union_mid_peak_overlapped -t $gsi_wo_reads -s hg38 -f bam -m -e 500 -o binding_csv/gsi_wo_ATAC_RPKM_peak_overlapped_union_e500.csv


####
python get_RPKM_on_regions_readcount.py -i $union_binding_mid -t $dmso_reads   -s hg38 -f bam -m -e 150 -o binding_csv/dmso_ATAC_RPKM_e150.csv
python get_RPKM_on_regions_readcount.py -i $union_binding_mid -t $gsi_3d_reads -s hg38 -f bam -m -e 150 -o binding_csv/gsi_3d_ATAC_RPKM_e150.csv
python get_RPKM_on_regions_readcount.py -i $union_binding_mid -t $gsi_wo_reads -s hg38 -f bam -m -e 150 -o binding_csv/gsi_wo_ATAC_RPKM_e150.csv

# python get_RPKM_on_regions_readcount.py -i $union_binding_mid -t $dmso_reads   -s hg38 -f bam -m -e 200 -o binding_csv/dmso_ATAC_RPKM_e200.csv
python get_RPKM_on_regions_readcount.py -i $union_binding_mid -t $gsi_3d_reads -s hg38 -f bam -m -e 200 -o binding_csv/gsi_3d_ATAC_RPKM_e200.csv
# python get_RPKM_on_regions_readcount.py -i $union_binding_mid -t $gsi_wo_reads -s hg38 -f bam -m -e 200 -o binding_csv/gsi_wo_ATAC_RPKM_e200.csv

# python get_RPKM_on_regions_readcount.py -i $union_binding_mid -t $dmso_reads   -s hg38 -f bam -m -e 500 -o binding_csv/dmso_ATAC_RPKM_e500.csv
python get_RPKM_on_regions_readcount.py -i $union_binding_mid -t $gsi_3d_reads -s hg38 -f bam -m -e 500 -o binding_csv/gsi_3d_ATAC_RPKM_e500.csv
# python get_RPKM_on_regions_readcount.py -i $union_binding_mid -t $gsi_wo_reads -s hg38 -f bam -m -e 500 -o binding_csv/gsi_wo_ATAC_RPKM_e500.csv

# python get_RPKM_on_regions_readcount.py -i $union_binding -t $dmso_reads   -s hg38 -f bam -o binding_csv/dmso_ATAC_RPKM_binding_site.csv
python get_RPKM_on_regions_readcount.py -i $union_binding -t $gsi_3d_reads -s hg38 -f bam -o binding_csv/gsi_3d_ATAC_RPKM_binding_site.csv
# python get_RPKM_on_regions_readcount.py -i $union_binding -t $gsi_wo_reads -s hg38 -f bam -o binding_csv/gsi_wo_ATAC_RPKM_binding_site.csv


