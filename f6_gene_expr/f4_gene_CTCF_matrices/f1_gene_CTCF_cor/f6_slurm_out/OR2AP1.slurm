#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o OR2AP1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i OR2AP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR2AT4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR2B11 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR2B2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR2B3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR2B6 
