#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=10000
#SBATCH -t 12:00:00
#SBATCH -p gpu 
#SBATCH -A zanglab
#SBATCH -o out_convert9.log

#Run program
python gene_CTCF_cor_ID_convert9.py
