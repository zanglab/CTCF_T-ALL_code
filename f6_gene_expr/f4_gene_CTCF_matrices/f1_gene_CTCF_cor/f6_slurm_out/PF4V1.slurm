#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PF4V1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PF4V1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PFAS 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PFDN1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PFDN2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PFDN4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PFDN5 
