#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o EFR3A.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i EFR3A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EFR3B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EFS 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EFTUD2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EGF 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EGFEM1P 
