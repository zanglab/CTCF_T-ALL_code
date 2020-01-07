#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ZNF793.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF793 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF799 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF80 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF800 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF804A 
