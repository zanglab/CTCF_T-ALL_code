#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ZNF165.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF165 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF169 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF17 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF174 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF175 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF177 
