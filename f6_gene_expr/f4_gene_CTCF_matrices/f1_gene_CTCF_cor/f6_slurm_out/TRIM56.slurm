#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TRIM56.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TRIM56 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TRIM58 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TRIM59 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TRIM6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TRIM6-TRIM34 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TRIM60 
