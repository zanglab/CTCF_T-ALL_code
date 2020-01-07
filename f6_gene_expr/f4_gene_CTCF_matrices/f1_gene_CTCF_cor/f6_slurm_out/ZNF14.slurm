#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ZNF14.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF14 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF140 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF141 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF142 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF143 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF146 
