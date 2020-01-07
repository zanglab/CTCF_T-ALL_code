#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ZNF764.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF764 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF765 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF766 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF767P 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF768 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF77 
