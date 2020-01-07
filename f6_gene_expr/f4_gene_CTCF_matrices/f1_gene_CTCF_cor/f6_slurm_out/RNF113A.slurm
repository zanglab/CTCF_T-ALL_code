#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o RNF113A.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i RNF113A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RNF113B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RNF114 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RNF115 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RNF121 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RNF122 
