#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ZNF18.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF18 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF180 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF181 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF182 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF184 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF185 
