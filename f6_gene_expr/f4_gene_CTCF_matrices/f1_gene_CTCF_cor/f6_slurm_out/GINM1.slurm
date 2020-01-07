#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o GINM1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i GINM1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GINS1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GINS2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GINS3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GINS4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GIP 
