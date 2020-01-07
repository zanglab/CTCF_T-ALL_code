#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ARL8A.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ARL8A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ARL8B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ARL9 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ARMC1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ARMC10 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ARMC12 
