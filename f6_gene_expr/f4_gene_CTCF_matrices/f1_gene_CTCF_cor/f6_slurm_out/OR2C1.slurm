#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o OR2C1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i OR2C1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR2C3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR2D2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR2D3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR2F1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR2F2 
