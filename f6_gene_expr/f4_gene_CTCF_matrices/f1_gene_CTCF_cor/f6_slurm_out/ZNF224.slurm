#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ZNF224.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF224 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF225 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF226 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF227 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF229 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF23 
