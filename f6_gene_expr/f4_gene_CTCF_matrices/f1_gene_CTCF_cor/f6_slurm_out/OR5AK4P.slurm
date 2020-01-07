#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o OR5AK4P.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i OR5AK4P 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR5AN1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR5AS1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR5AU1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR5B12 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR5B2 
