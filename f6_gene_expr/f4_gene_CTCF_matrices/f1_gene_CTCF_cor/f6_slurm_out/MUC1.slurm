#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o MUC1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i MUC1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MUC12 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MUC13 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MUC15 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MUC16 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MUC17 
