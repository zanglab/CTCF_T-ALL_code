#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o RNF25.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i RNF25 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RNF26 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RNF31 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RNF32 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RNF34 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RNF38 
