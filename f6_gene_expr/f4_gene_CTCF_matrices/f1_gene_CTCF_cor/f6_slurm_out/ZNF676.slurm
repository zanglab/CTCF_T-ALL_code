#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ZNF676.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF676 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF677 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF678 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF679 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF680 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF681 
