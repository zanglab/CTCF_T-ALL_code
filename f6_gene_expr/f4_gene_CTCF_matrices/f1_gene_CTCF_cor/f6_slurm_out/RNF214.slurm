#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o RNF214.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i RNF214 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RNF215 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RNF216 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RNF216P1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RNF217 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RNF219 
