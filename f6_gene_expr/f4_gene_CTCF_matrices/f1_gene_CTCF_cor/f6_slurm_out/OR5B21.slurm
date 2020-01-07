#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o OR5B21.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i OR5B21 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR5C1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR5E1P 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR5H1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR5H14 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR5H6 
