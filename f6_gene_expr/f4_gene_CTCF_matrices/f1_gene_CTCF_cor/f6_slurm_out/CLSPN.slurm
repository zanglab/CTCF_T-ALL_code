#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CLSPN.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CLSPN 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CLSTN1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CLSTN2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CLSTN3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CLTA 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CLTB 
