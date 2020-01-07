#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o RAD18.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i RAD18 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RAD21 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RAD21L1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RAD23A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RAD23B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RAD50 
