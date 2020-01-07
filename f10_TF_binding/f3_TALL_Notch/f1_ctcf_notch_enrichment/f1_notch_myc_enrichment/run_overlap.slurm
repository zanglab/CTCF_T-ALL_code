#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=40000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o py1.out
#Run program


domain_bed="/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f6_gene_expr/f4_gene_CTCF_matrices/f2_gene_CTCF_domain/all_CTCF_domainInfo_GT100K_LT1M_EachSide.bed"


notch_file_dir="/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/12_new_gained_lost_from_all_features/f8_cofactor_T-ALL_Notch/T_ALL/f0_dynamic_notch_co_binding/Cutll1_Notch1_files" 
time python find_overlap_keep_info_NOT_sep_strand_lastColMarked_notch_revised.py -a $domain_bed -b ${notch_file_dir}/GSM1252936_dynamic_NOTCH1_peaks_hg38.bed -s hg38 -p dynamic_notch_domain_co_occurrence.csv
time python find_overlap_keep_info_NOT_sep_strand_lastColMarked.py -a $domain_bed -b ${notch_file_dir}/CUTLL1_NOTCH1_all_summits.bed -s hg38 -p cutll1_notch_domain_co_occurrence.csv


myc_file_dir="/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/18_MYC-ChIP_shCTCF-RNA/MYC_ChIP/f0_processing/callpeak/macs2_out"
time python find_overlap_keep_info_NOT_sep_strand_lastColMarked.py -a $domain_bed -b ${myc_file_dir}/CUTLL1_MYC_summits.bed -s hg38 -p cutll1_MYC_domain_co_occurrence.csv
