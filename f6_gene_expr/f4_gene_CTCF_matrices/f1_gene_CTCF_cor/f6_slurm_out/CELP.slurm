#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CELP.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CELP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CELSR1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CELSR2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CELSR3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CEMIP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CEMP1 
