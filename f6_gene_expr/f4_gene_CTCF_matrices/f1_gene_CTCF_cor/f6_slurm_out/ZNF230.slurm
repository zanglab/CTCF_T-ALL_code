#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ZNF230.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF230 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF232 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF233 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF234 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF235 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF236 
