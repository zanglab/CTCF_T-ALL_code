#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ZNF492.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF492 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF493 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF496 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF497 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF500 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF501 
