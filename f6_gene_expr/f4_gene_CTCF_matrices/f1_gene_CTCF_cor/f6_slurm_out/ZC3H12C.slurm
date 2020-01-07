#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ZC3H12C.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ZC3H12C 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZC3H12D 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZC3H13 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZC3H14 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZC3H15 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZC3H18 
