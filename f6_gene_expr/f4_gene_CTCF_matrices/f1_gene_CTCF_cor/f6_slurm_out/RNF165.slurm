#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o RNF165.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i RNF165 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RNF166 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RNF167 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RNF168 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RNF169 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RNF17 
