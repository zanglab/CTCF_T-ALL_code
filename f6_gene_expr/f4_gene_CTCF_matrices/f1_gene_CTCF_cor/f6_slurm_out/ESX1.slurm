#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ESX1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ESX1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ESYT1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ESYT2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ESYT3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ETAA1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ETF1 
