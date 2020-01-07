#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o MAP7D1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i MAP7D1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MAP7D2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MAP7D3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MAP9 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MAPK1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MAPK10 
