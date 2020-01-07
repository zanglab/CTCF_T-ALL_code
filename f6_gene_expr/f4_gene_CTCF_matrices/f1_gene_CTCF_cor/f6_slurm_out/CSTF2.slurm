#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CSTF2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CSTF2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CSTF2T 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CSTF3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CSTL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CT45A1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CT45A10 
