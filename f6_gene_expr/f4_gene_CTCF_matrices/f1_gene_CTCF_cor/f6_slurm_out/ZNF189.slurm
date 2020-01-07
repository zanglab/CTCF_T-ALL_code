#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ZNF189.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF189 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF19 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF192P1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF195 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF197 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF2 
