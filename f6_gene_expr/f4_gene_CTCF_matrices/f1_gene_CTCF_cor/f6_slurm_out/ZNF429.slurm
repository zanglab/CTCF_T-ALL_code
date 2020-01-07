#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ZNF429.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF429 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF43 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF430 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF431 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF432 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF433 
