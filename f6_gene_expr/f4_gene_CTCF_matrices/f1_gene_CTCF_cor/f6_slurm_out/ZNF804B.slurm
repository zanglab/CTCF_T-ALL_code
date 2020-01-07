#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ZNF804B.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF804B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF805 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF806 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF808 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF81 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF813 
