#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o OR6K2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i OR6K2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR6K3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR6M1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR6N1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR6N2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR6P1 
