#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PCDHA11.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PCDHA11 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PCDHA12 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PCDHA13 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PCDHA2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PCDHA3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PCDHA4 
