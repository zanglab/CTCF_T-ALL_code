#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o OR2L5.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i OR2L5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR2L8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR2M1P 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR2M3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR2M4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR2S2 
