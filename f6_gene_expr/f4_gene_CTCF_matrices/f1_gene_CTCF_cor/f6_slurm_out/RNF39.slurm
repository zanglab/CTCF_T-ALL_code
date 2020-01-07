#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o RNF39.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i RNF39 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RNF4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RNF40 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RNF41 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RNF43 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RNF44 
