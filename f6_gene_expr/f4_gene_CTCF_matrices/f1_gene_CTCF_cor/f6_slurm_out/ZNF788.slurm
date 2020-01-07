#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ZNF788.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF788 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF789 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF79 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF790 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF791 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF792 
