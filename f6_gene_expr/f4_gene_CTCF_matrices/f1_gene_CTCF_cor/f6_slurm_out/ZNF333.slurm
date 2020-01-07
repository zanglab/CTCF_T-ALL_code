#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ZNF333.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF333 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF334 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF335 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF337 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF33A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF33B 
