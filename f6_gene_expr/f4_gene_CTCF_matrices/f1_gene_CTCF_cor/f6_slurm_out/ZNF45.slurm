#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ZNF45.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF45 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF451 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF454 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF460 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF461 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF462 
