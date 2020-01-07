#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ITGB4.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ITGB4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ITGB5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ITGB6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ITGB7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ITGB8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ITGBL1 
