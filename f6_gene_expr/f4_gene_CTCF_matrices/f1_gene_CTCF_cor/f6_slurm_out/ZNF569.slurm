#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ZNF569.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF569 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF57 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF570 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF571 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF572 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF573 
