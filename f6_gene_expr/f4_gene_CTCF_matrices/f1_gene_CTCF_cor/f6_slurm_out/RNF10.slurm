#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o RNF10.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i RNF10 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RNF103 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RNF103-CHMP3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RNF11 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RNF111 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RNF112 
