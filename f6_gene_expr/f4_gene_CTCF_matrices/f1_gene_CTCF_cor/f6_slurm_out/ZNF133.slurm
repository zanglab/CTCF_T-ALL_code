#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ZNF133.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF133 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF134 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF135 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF136 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF137P 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF138 
