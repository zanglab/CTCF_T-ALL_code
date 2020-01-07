#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ZNF347.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF347 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF35 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF350 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF354A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF354B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF354C 
