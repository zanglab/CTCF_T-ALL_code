#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o GJA1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i GJA1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GJA10 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GJA3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GJA4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GJA5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GJA8 
