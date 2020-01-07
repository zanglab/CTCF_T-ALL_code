#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TRPM4.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TRPM4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TRPM5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TRPM6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TRPM7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TRPM8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TRPS1 
