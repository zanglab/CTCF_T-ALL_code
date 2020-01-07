#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o RNF20.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i RNF20 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RNF207 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RNF208 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RNF212 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RNF212B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RNF213 
