#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ARL3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ARL3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ARL4A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ARL4C 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ARL4D 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ARL5A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ARL5B 
