

=======
Notes
=======

- patient_label = {'PD9':'Patient1', 'PD31 or PD30':'Patient2',}
- For PRAD, results are using the one comparing PRAD vs. 4 prostate tissues (PRAD_TissueAdded)
- For TCGA related analysis (chromatin accessibility and clinical outcome), TCGA Pan-cancer peaks are used 
- Max of Days to last follow up and Days to death was used in survival analysis.
- Everything keep to CUTLL1 (NOTCH ChIP, gene expression, Hi-C) for shCTCF associated results


====
union CTCF
====

- occupancy score vs. % motif 
	- dir	./fz_figures/union_binding_features/unionCTCF_occupancy_score_vs_motif.py

- occupancy score power model fit 
	- dir	./fz_figures/union_binding_features/unionCTCF_occupancy_score_power_model_fit.py

- distribution of number of CTCF peaks 
	- dir	./fz_figures/union_binding_features/unionCTCF_peak_distribution.py

- distribution of intervals of union CTCF peaks
	- dir	./fz_figures/union_binding_features/unionCTCF_summit_intervals.py

- union peak summit intervals
	- dir	./f2_union_CTCFs/f1_union_all_summits/py5_plot_summit_intervals.py
	
- pan Cancer CTCF overlap
	- dir	./f3_cancer_specific_CTCFs/f6_panCancer_specific_CTCF_overlap
	
- pan Cancer CTCF annotation/location
	- dir	./f3_cancer_specific_CTCFs/f7_panCancer_specific_CTCF_location
	
- CTCF binding signal of lost vs. constitutive vs. gained in T-ALL
	- dir	./f3_cancer_specific_CTCFs/f5_binding_signal_compr/f2_figs


====
CTCF heatmap
====

- CTCF heatmap of T-ALL cell lines with SMC3 binding
	- dir	./fz_figures/ctcf_TALL_with_SMC3

- CTCF heatmap T-ALL patients
	- dir ./f3_cancer_specific_CTCFs/f8_panCancer_CTCF_binding_pattern/f3_cancerType_CTCF_heatmap/
	
- CTCF heatmap of Pan-cancers with DNA me. change 
	- dir	./fz_figures/ctcf_heatmap_DNAme_mutation/pancancer_fig_DNAme

- CTCF heatmap of Pan-cancers with with mutation
	dir	./fz_figures/ctcf_heatmap_DNAme_mutation/{pancancer_fig_mutation


====
ATAC-seq
====

- TCGA pan Cancer ATAC
	- dir	./f5_TCGA_ATAC/f2_cancer_specific_CTCF_accessibility_panCancer/f3_sig_compr_layout_fc/mynorm
	
- T-ALL ATAC
	- dir	./f5_TCGA_ATAC/f6_TALL_ATAC_raw_count
	

====
Survival
====

- TCGA survival
	- dir	./f5_TCGA_ATAC/f4_clinical_survival_panCancer
	

====
Hi-C
====
	
- volcano plot & box plot
	- dir	./f4_HiC_TADs/f2_gained_lost_hic_interaction/f4_scatter_with_boxplot_figs
	
- mirror bar plot
	- dir	./f4_HiC_TADs/f2_gained_lost_hic_interaction/f5_mirror_bar_figs


====
histone modification
====

- pan cancer
	- dir	./f9_histone_modification/f3_HM_heatmap/f3_HM_heatmap_combined_PanCancer

- T-ALL patient samples
	- dir	./f9_histone_modification/f3_HM_heatmap/f4_HM_heatmap_combined_patients


====
mutation
====

- pan cancer (bar plot centered at CTCF +/-200bp)
	- dir	./f7_ICGC_mutation/panCancer_ICGC/f5_from_motif_mid_pos_replot_f3_figs
	
- T-ALL (bar plot centered at CTCF)
	- dir	./f7_ICGC_mutation/T_ALL_new/f5_from_motif_mid_pos_replot_f3_figs


====
Hi-C TADs compared with constitutive CTCF
====

- number of CTCF & number of merged region
	- dir	./f4_HiC_TADs/f4_conserved_CTCF_loops_compr_TADs/*.bed
-  number of overlapped
	- dir	./f4_HiC_TADs/f4_conserved_CTCF_loops_compr_TADs/overlapped/*merged*bed


====
correlation of CTCF-gene pairs
====

- scatter plot of CTCF binding vs. gene expr
	- dir	./f6_gene_expr/f5_gene_CTCF_cor_figs/ctcf_gene_cor_scatter/

- genome wide CTCF gene inter/intra-domain correlation
	- dir	./f6_gene_expr/f5_gene_CTCF_cor_figs/genome_intra_inter_domain_cor_compr
	
- pan cancer CTCF intra-domain hi-cor genes
	- dir	./f6_gene_expr/f7_panCancer_deg_figs/f7_panCancer_intra_domain_hicor_figs/

- T-ALL CTCF intra-domain gene expression patterns, by DNA methylation
	- dir	./f6_gene_expr/f7_panCancer_deg_figs/f5_TALL_patient_intra_domain_deg_pattern_figs/

- pan Cancer CTCF intra-domain gene expression patterns
	- dir	./f6_gene_expr/f7_panCancer_deg_figs/f4_panCancer_intra_domain_deg_pattern_figs/

- pan cancer CTCF intra-domain gene expression patterns, by DNA methylation
	- dir	./f6_gene_expr/f7_panCancer_deg_figs/f4_panCancer_intra_domain_deg_pattern_figs_by_DNAme/
	
	
====
T-ALLgained NOTCH1 enrichment
====

- BART prediction (intra-domain of GT100K and LT1M)
	- dir	./f10_TF_binding/f2_panCancer_intra_domain_GT100K_LT1M_TF_enrichment/f4_scatter_plot/

- CTCF~NOTCH1 Hi-C interaction z-scores
	- dir	./f10_TF_binding/f3_TALL_Notch/f6_CTCF_dynamic_notch_targetGene/f4_interaction_compr_zscores

- enrichment of NOTCH1 
	- dir	./f10_TF_binding/f3_TALL_Notch/f1_ctcf_notch_enrichment/f2_figs_notch_enrichment_layouts

- enrichment of NOTCH1 and MYC
	- dir	./f10_TF_binding/f3_TALL_Notch/f1_ctcf_notch_enrichment/f3_notch_myc_enrichment_figs

- gained CTCF ~ H3K27ac
	- dir	./f10_TF_binding/f3_TALL_Notch/f2_notch_H3K27ac

- CTCF ~ NOTCH1 pairwise distances
	- dir	./f10_TF_binding/f3_TALL_Notch/f4_ctcf_notch_pairwise_dis/f2_dis_figs


====
NOTCH1, BAF, CTCF
====

- CTCF binding changes in GSI treatment
	- dir	./f11_GSI_shCTCF/f1_GSI_chip_binding/f3_tall_gained_GSI_bidings_changes

- compare ATAC-seq RPKM of T-ALL gained CTCFs in DMSO, GSI and GSI washout
	- dir	./f11_GSI_shCTCF/f5_GSI_ATAC_chip_20191029/f1_results_from_bam/f3_gained_CTCF_RPKM_compr/

- shCTCF DEG (down genes)  vs. in T-ALL gained CTCFs targets (log2FC<-0.26, FDR<0.001, #987)
	- dir	./f11_GSI_shCTCF/f2_shCTCF_TALL_gained/f2_figs

- BART prediction of shCTCF DEG (down genes) (log2FC<-0.58, FDR<0.01, #674)
	- dir	./f11_GSI_shCTCF/shCTCF_RNA/f3_shCTCF_DEG_BART/rank_dot_figs/

- shCTCF DEG (down genes) vs. NOTCH1 target genes
	dir	./f11_GSI_shCTCF/shCTCF_RNA/f5_shCTCF_DEG_vs_NOTCH1_targets_vs_TALL_DEG/
	
- shCTCF DEG (down genes) vs. GSI-down, GSI-wo-up genes
	- dir	./f11_GSI_shCTCF/f6_Jurkat_GSI_RNA_seq_new/f5_enrichment_of_dynamic_genes/

- expression in CUTLL1 vs. T-cell of dyNOTCH1 (and gained CTCF) intra-domain genes
	- dir	./f10_TF_binding/f3_TALL_Notch/f9_dynamicNOTCH_gained_CTCF_gene_expr/f2_dNOTCH_gainedCTCF_figs

- expression of CUTLL1 up-regulated NOTCH1+gained CTCF intra-domain genes in shCTCF, by dyNOTCH ~ gained CTCF Hi-C interaction
	- dir	./f10_TF_binding/f3_TALL_Notch/f9_dynamicNOTCH_gained_CTCF_gene_expr/f4_dNOTCH_gainedCTCF_expr_in_shCTCF_figs

- SWI/SNF, T-ALL BRG1/CTCF on T-ALL gained
	- dir	./f12_SWI_SNF/f3_CUTLL1_BRG1/f1_brg1_binding_compr/f2_chip_binding_compr_CTCF_types
	
- AML SMARCA4/CTCF binding on AML gained
	- dir	./f12_SWI_SNF/f2_AML_site_plot/f2_chip_binding_compr_CTCF_types
	
- BRCA/PRAD Crispr screens
	- dir	./fz_figures/ctcf_crispr_screens_PMID31727847


