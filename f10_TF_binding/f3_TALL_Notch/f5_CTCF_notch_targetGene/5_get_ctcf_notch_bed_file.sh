mkdir f5_ctcf_notch_bed_file
tail -n +2 f4_interaction_compr_zscores/raw_T-ALL_gained_zscore_delta_change_filtered.txt|awk '{OFS="\t";print$2,$3-10,$3+10,$1}'> f5_ctcf_notch_bed_file/gaiend_CTCF.be
tail -n +2 f4_interaction_compr_zscores/raw_T-ALL_gained_zscore_delta_change_filtered.txt|awk '{OFS="\t";print$2,$4-10,$4+10,"NOTCH",$NF}'> f5_ctcf_notch_bed_file/gaiend_CTCF_hic_increased_notch.bed
