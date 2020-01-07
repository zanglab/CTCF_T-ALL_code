
atac_dir="/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/12_new_gained_lost_from_all_features/f13_ATAC_Jurkat_CD4/f1_proprocess/"

cat ${atac_dir}/GSM2411156/GSM2411156_Jurkat_Untreated_R1_ATAC-seq_peaks.narrowPeak \
${atac_dir}/GSM2411157/GSM2411157_Jurkat_Untreated_R2_ATAC-seq_peaks.narrowPeak \
${atac_dir}/GSM2411158/GSM2411158_Jurkat_Untreated_R3_ATAC-seq_peaks.narrowPeak > Jurkat_ATAC_peaks.narrowPeak

cat ${atac_dir}/GSM2411156/GSM2411156_Jurkat_Untreated_R1_ATAC-seq_PEnoM.bed \
${atac_dir}/GSM2411157/GSM2411157_Jurkat_Untreated_R2_ATAC-seq_PEnoM.bed \
${atac_dir}/GSM2411158/GSM2411158_Jurkat_Untreated_R3_ATAC-seq_PEnoM.bed  > Jurkat_ATAC_PEnoM.bed


python find_overlap_keep_info_NOT_sep_strand.py -a data/union_CTCF_No_Jurkat_DNAmethylation_with_motif.bed -b data/Jurkat_ATAC_peaks.narrowPeak -p data/union_CTCF_No_Jurkat_DNAmethylation_with_motif_ATAC_peak_overlapped.bed -s hg38
# 
# cat union_CTCF_No_Jurkat_DNAmethylation_with_motif_ATAC_peak_overlapped.bed|awk '{OFS="\t";print$1,$NF,$NF+1,$4}' > union_CTCF_No_Jurkat_DNAmethylation_with_motif_ATAC_peak_overlapped_mid_position.bed


