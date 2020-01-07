#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ZNF20.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF20 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF200 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF202 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF204P 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF205 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF207 
