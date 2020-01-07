#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o RNF5.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i RNF5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RNF5P1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RNF6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RNF7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RNF8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RNFT1 
