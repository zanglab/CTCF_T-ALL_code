#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SMARCA4.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SMARCA4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SMARCA5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SMARCAD1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SMARCAL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SMARCB1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SMARCC1 
