#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=60000
#SBATCH -t 24:00:00
#SBATCH -p largemem
#SBATCH -A zanglab

#Run program



infile="data/ctcf_cor_atac_union_site_NoDNAme_withMotif.bed"
bart2 region -i ${infile} -s hg38 --outdir bart2_out

infile="data/ctcf_cor_atac_union_site_NoDNAme_withMotif_JurkatCTCF.bed"
bart2 region -i ${infile} -s hg38 --outdir bart2_out

infile="data/ctcf_cor_atac_union_site_NoDNAme_withMotif_NoJurkatCTCF.bed"
bart2 region -i ${infile} -s hg38 --outdir bart2_out

infile="data/union_CTCF_No_Jurkat_DNAmethylation_with_motif_ATAC_peak_overlapped.bed"
bart2 region -i ${infile} -s hg38 --outdir bart2_out

infile="data/union_CTCF_No_Jurkat_DNAmethylation_with_motif_ATAC_peak_overlapped_Jurkat_CTCF_NOT_overlapped.bed"
bart2 region -i ${infile} -s hg38 --outdir bart2_out

infile="data/union_CTCF_No_Jurkat_DNAmethylation_with_motif_ATAC_peak_overlapped_Jurkat_CTCF_overlapped.bed"
bart2 region -i ${infile} -s hg38 --outdir bart2_out


