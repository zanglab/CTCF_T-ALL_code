#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ZNF208.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF208 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF211 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF212 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF213 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF214 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF215 
