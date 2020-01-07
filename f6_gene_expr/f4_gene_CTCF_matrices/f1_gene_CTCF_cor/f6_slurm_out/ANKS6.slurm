#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ANKS6.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ANKS6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ANKUB1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ANKZF1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ANLN 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ANO1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ANO10 
