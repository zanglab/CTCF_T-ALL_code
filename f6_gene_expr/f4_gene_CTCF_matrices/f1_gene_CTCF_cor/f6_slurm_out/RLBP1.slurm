#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o RLBP1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i RLBP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RLF 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RLIM 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RLN1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RLN2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RLN3 
