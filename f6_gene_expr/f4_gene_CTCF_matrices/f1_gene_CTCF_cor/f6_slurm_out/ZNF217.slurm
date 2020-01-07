#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ZNF217.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF217 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF219 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF22 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF221 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF222 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF223 
