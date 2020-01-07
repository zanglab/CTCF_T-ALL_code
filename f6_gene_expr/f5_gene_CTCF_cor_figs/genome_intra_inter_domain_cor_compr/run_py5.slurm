#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=60000
#SBATCH -t 24:00:00
#SBATCH -p largemem 
#SBATCH -A zanglab
#SBATCH -o py5.out
#Run program
time python 5_intra_inter_gene_CTCF_cor_domain_GT100k_LT1M.py

