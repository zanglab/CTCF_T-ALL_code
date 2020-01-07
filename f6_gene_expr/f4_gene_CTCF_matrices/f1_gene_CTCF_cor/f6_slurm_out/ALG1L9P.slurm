#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ALG1L9P.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ALG1L9P 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ALG2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ALG3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ALG5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ALG6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ALG8 
