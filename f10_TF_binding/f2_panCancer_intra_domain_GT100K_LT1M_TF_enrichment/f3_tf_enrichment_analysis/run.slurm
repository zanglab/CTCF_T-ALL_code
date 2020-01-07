#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=20000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o slurm_out.log



# python run_udhs_motif_enrichment_hg38.py -t data/T-ALL_gained.bed  -o test1_out
# time python run_udhs_motif_enrichment_hg38.py -t data/T-ALL_gained_raw_center_hicor_promoter_overlapped.bed -c data/T-ALL_gained_raw_all_ctrl_promoter_overlapped.bed -o test2_out 

python run_udhs_motif_enrichment_hg38.py -t ../f2_high_cor_high_logFC_increased_decreased/T-ALL_gained_center_hicor.bed -c ../f2_high_cor_high_logFC_increased_decreased/T-ALL_gained_center_ctrl.bed -o tf_enrichment_results
python run_udhs_motif_enrichment_hg38.py -t ../f2_high_cor_high_logFC_increased_decreased/T-ALL_lost_center_hicor.bed -c ../f2_high_cor_high_logFC_increased_decreased/T-ALL_lost_center_ctrl.bed -o tf_enrichment_results
python run_udhs_motif_enrichment_hg38.py -t ../f2_high_cor_high_logFC_increased_decreased/CRC_gained_center_hicor.bed -c ../f2_high_cor_high_logFC_increased_decreased/CRC_gained_center_ctrl.bed -o tf_enrichment_results
python run_udhs_motif_enrichment_hg38.py -t ../f2_high_cor_high_logFC_increased_decreased/CRC_lost_center_hicor.bed -c ../f2_high_cor_high_logFC_increased_decreased/CRC_lost_center_ctrl.bed -o tf_enrichment_results
