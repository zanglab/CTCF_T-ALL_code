#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o OR2T29.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i OR2T29 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR2T3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR2T33 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR2T34 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR2T35 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR2T4 
