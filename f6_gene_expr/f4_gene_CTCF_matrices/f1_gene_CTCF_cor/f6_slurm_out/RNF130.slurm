#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o RNF130.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i RNF130 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RNF133 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RNF135 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RNF138 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RNF138P1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RNF139 
