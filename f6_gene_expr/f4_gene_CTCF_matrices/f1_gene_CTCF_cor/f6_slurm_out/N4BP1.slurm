#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o N4BP1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i N4BP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i N4BP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i N4BP2L1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i N4BP2L2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i N4BP3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i N6AMT1 
