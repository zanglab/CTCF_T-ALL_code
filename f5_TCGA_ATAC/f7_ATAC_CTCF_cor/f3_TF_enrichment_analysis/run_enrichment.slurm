#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=40000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab

#Run program

# python run_udhs_motif_enrichment_hg38.py -t data/T-ALL_gained.bed  -o test1_out

# time python run_udhs_motif_enrichment_hg38.py -t data/ctcf_cor_atac_union_site_NoDNAme_withMotif.bed -o TF_enrichemnt
# time python run_udhs_motif_enrichment_hg38.py -t data/ctcf_cor_atac_union_site_NoDNAme_withMotif_JurkatCTCF.bed -o TF_enrichemnt
# time python run_udhs_motif_enrichment_hg38.py -t data/ctcf_cor_atac_union_site_NoDNAme_withMotif_NoJurkatCTCF.bed -o TF_enrichemnt
# time python run_udhs_motif_enrichment_hg38.py -t data/union_CTCF_No_Jurkat_DNAmethylation_with_motif_ATAC_peak_overlapped.bed -o TF_enrichemnt
# time python run_udhs_motif_enrichment_hg38.py -t data/union_CTCF_No_Jurkat_DNAmethylation_with_motif_ATAC_peak_overlapped_Jurkat_CTCF_NOT_overlapped.bed -o TF_enrichemnt
# time python run_udhs_motif_enrichment_hg38.py -t data/union_CTCF_No_Jurkat_DNAmethylation_with_motif_ATAC_peak_overlapped_Jurkat_CTCF_overlapped.bed -o TF_enrichemnt

t_file="data/union_CTCF_No_Jurkat_DNAmethylation_with_motif_ATAC_peak_overlapped_Jurkat_CTCF_NOT_overlapped.bed"
c_file="data/union_CTCF_No_Jurkat_DNAmethylation_with_motif_ATAC_peak_overlapped_Jurkat_CTCF_overlapped.bed"
time python run_udhs_motif_enrichment_hg38.py -t ${t_file} -c ${c_file} -o TF_enrichemnt_with_control

t_file="data/ctcf_cor_atac_union_site_NoDNAme_withMotif_NoJurkatCTCF.bed"
c_file="data/ctcf_cor_atac_union_site_NoDNAme_withMotif_JurkatCTCF.bed"
time python run_udhs_motif_enrichment_hg38.py -t ${t_file} -c ${c_file} -o TF_enrichemnt_with_control
