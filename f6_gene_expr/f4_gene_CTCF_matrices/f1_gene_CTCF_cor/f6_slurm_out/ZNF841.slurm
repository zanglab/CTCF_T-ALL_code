#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ZNF841.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF841 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF843 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF844 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF845 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF846 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF847P 
