#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ZNF252P.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF252P 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF253 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF254 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF256 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF257 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF26 
