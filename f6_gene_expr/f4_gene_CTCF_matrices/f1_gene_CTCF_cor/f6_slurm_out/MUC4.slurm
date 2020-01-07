#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o MUC4.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i MUC4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MUC5AC 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MUC5B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MUC6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MUC7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MUCL1 
