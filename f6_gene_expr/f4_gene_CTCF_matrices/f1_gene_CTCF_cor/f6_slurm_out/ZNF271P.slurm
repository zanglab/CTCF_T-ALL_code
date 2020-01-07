#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ZNF271P.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF271P 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF273 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF274 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF275 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF276 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF277 
