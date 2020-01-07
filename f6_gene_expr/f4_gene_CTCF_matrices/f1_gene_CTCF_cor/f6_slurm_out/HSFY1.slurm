#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o HSFY1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i HSFY1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HSFY1P1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HSFY2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HSH2D 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HSP90AA1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HSP90AB1 
