#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ZNF580.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF580 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF581 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF582 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF583 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF584 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF585A 
