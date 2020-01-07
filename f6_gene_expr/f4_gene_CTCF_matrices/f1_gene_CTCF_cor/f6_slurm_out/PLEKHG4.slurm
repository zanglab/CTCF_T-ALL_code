#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PLEKHG4.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PLEKHG4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PLEKHG4B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PLEKHG5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PLEKHG6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PLEKHG7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PLEKHH1 
