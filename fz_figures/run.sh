## ==== union CTCFs 
python unionCTCF_peak_distribution.py
python unionCTCF_summit_intervals.py 
python unionCTCF_occupancy_score_power_model_fit.py > log/unionCTCF_occupancy_score_power_model_fit.log
python unionCTCF_occupancy_score_vs_motif.py




## ==== cancer_specific
python cancer_specific_occupancy_heatmap.py 
python cancer_specific_binding_signal_compr_TALL_gained_lost_constitutive.py
python cancer_specific_binding_signal_compr_panCancer_gained_constitutive.py
python cancer_specific_binding_overlap.py
python cancer_specific_binding_gene_annotation.py 

## ==== TCGA-ATAC
python cancer_specific_TCGA_ATAC_accessibility_compr.py 


## ==== gene expression ~ CTCF binding
python gene_expr_CTCF_intra_inter_hicor_compr.py > log/gene_expr_CTCF_intra_inter_hicor_compr.log


