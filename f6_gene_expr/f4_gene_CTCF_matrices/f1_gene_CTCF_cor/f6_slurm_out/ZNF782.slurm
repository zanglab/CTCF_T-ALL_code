#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ZNF782.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF782 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF783 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF784 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF785 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF786 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF787 
