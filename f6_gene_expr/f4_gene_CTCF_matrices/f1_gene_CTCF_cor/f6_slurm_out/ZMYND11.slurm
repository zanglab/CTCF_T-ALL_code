#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ZMYND11.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ZMYND11 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZMYND12 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZMYND15 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZMYND19 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZMYND8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF10 
