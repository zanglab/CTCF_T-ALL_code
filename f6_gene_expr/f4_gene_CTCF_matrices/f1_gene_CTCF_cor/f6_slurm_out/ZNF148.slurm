#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ZNF148.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF148 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF154 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF155 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF157 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF16 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF160 
