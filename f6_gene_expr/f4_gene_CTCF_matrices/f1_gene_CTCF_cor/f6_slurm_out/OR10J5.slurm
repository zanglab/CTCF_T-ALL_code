#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o OR10J5.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i OR10J5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR10K1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR10K2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR10P1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR10Q1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR10R2 
