#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o POTEKP.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i POTEKP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i POTEM 
time python ../7_cor_CTCF_binding_GeneExpr.py -i POU1F1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i POU2AF1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i POU2F1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i POU2F2 
