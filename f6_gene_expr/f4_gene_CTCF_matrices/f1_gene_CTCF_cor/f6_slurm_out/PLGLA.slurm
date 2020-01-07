#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PLGLA.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PLGLA 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PLGLB1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PLGLB2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PLGRKT 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PLIN1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PLIN2 
