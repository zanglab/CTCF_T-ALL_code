#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ZNF526.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF526 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF527 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF528 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF529 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF530 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF532 
