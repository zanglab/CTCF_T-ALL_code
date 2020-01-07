#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=40000
#SBATCH -t 24:00:00
#SBATCH -p gpu
#SBATCH -A zanglab
#SBATCH -o out_overlap.log

cancerType="BRCA"
python find_overlap_keep_info_NOT_sep_strand_revised.py -a mynorm_TCGA-ATAC_PanCan_Log2_QuantileNorm_Counts_plus5.txt -b TCGA-ATAC_Cancer_Type-specific_PeakCalls/${cancerType}_peakCalls.txt -s hg38 -p overlapped/ATAC_PanCan_Log2_QN_${cancerType}_overlapped.txt -q overlapped/ATAC_PanCan_Log2_QN_${cancerType}_NONoverlapped.txt

cancerType="COAD"
python find_overlap_keep_info_NOT_sep_strand_revised.py -a mynorm_TCGA-ATAC_PanCan_Log2_QuantileNorm_Counts_plus5.txt -b TCGA-ATAC_Cancer_Type-specific_PeakCalls/${cancerType}_peakCalls.txt -s hg38 -p overlapped/ATAC_PanCan_Log2_QN_${cancerType}_overlapped.txt -q overlapped/ATAC_PanCan_Log2_QN_${cancerType}_NONoverlapped.txt

cancerType="LUAD"
python find_overlap_keep_info_NOT_sep_strand_revised.py -a mynorm_TCGA-ATAC_PanCan_Log2_QuantileNorm_Counts_plus5.txt -b TCGA-ATAC_Cancer_Type-specific_PeakCalls/${cancerType}_peakCalls.txt -s hg38 -p overlapped/ATAC_PanCan_Log2_QN_${cancerType}_overlapped.txt -q overlapped/ATAC_PanCan_Log2_QN_${cancerType}_NONoverlapped.txt

cancerType="PRAD"
python find_overlap_keep_info_NOT_sep_strand_revised.py -a mynorm_TCGA-ATAC_PanCan_Log2_QuantileNorm_Counts_plus5.txt -b TCGA-ATAC_Cancer_Type-specific_PeakCalls/${cancerType}_peakCalls.txt -s hg38 -p overlapped/ATAC_PanCan_Log2_QN_${cancerType}_overlapped.txt -q overlapped/ATAC_PanCan_Log2_QN_${cancerType}_NONoverlapped.txt

