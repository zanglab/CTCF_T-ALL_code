#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PMEL.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PMEL 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PMEPA1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PMF1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PMF1-BGLAP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PMFBP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PML 
