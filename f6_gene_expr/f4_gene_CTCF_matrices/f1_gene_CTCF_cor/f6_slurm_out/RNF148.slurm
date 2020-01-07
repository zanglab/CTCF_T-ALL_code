#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o RNF148.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i RNF148 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RNF149 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RNF150 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RNF151 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RNF152 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RNF157 
