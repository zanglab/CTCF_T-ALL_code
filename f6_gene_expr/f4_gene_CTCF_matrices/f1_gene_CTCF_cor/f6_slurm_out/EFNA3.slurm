#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o EFNA3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i EFNA3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EFNA4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EFNA5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EFNB1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EFNB2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EFNB3 
