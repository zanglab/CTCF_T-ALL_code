#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ZNF397.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF397 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF398 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF404 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF407 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF408 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF41 
