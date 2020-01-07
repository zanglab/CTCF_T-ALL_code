#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o MYL2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i MYL2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MYL3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MYL4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MYL5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MYL6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MYL6B 
