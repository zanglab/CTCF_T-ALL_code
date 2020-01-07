#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ZNF33BP1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF33BP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF34 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF341 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF343 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF345 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF346 
