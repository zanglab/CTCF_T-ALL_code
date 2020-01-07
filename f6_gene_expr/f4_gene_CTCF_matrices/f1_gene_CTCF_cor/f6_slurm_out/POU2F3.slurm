#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o POU2F3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i POU2F3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i POU3F1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i POU3F2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i POU3F3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i POU3F4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i POU4F1 
