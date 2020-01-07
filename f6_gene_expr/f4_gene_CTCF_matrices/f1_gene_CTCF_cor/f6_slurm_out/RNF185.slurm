#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o RNF185.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i RNF185 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RNF186 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RNF187 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RNF19A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RNF19B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RNF2 
