#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o OR2G2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i OR2G2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR2G3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR2G6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR2H1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR2H2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR2J2 
