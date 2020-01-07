#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o OR1N2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i OR1N2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR1Q1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR1S1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR2A1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR2A12 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR2A14 
