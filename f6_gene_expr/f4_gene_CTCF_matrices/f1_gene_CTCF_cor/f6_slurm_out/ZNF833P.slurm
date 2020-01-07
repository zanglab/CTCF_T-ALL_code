#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ZNF833P.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF833P 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF835 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF836 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF837 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF839 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF84 
