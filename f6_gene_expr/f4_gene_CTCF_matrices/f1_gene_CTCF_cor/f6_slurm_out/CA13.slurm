#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CA13.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CA13 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CA14 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CA2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CA3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CA4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CA5A 
