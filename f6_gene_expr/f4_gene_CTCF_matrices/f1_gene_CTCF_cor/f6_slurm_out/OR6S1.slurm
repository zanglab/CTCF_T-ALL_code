#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o OR6S1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i OR6S1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR6T1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR6V1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR6W1P 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR6X1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR6Y1 
