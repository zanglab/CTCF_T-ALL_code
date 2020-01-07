
ctcf_dir="/nv/vol190/zanglab/zw5j/projects_data/T_ALL_CTCF_panos/CTCF/peaks"

cat ${ctcf_dir}/JURKAT_CTCF1_peaks.narrowPeak \
${ctcf_dir}/JURKAT_CTCF2_peaks.narrowPeak > data/Jurkat_CTCF_peaks.narrowPeak


a_file="data/union_CTCF_No_Jurkat_DNAmethylation_with_motif_ATAC_peak_overlapped.bed"
b_file="data/Jurkat_CTCF_peaks.narrowPeak"

python find_overlap_keep_info_NOT_sep_strand.py -a ${a_file} -b ${b_file} \
-p data/union_CTCF_No_Jurkat_DNAmethylation_with_motif_ATAC_peak_overlapped_Jurkat_CTCF_overlapped.bed \
-q data/union_CTCF_No_Jurkat_DNAmethylation_with_motif_ATAC_peak_overlapped_Jurkat_CTCF_NOT_overlapped.bed -s hg38
# 
# cat union_CTCF_No_Jurkat_DNAmethylation_with_motif_ATAC_peak_overlapped.bed|awk '{OFS="\t";print$1,$NF,$NF+1,$4}' > union_CTCF_No_Jurkat_DNAmethylation_with_motif_ATAC_peak_overlapped_mid_position.bed


