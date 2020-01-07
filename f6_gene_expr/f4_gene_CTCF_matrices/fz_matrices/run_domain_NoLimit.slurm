#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=40000
#SBATCH -t 24:00:00
#SBATCH -p standard 
#SBATCH -A zanglab
#SBATCH -o slurm_domain_NoLimit.out
time python return_gene_CTCF_same_domain_NoLimit_prediction.py
