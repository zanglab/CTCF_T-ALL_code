#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ZNF599.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF599 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF600 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF605 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF606 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF607 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF608 
