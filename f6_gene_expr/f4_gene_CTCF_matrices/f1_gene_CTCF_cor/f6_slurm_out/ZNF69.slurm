#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ZNF69.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF69 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF691 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF692 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF695 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF696 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF697 
