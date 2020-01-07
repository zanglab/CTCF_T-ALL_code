#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ZC3H3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ZC3H3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZC3H4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZC3H6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZC3H7A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZC3H7B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZC3H8 
