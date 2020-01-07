#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ITGB1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ITGB1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ITGB1BP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ITGB1BP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ITGB2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ITGB3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ITGB3BP 
