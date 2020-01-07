#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=40000
#SBATCH -t 24:00:00
#SBATCH -p largemem 
#SBATCH -A zanglab
#SBATCH -o slurm_domain_100k_1M.out
time python return_gene_CTCF_same_domain_GT100k_LT1M_prediction.py
