#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ZNF826P.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF826P 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF827 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF829 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF83 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF830 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF831 
