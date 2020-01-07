#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o RNF14.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i RNF14 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RNF141 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RNF144A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RNF144B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RNF145 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RNF146 
