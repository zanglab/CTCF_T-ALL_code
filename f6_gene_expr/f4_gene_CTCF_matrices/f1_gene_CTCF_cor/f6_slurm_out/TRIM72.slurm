#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TRIM72.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TRIM72 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TRIM73 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TRIM74 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TRIM77 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TRIM8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TRIM9 
