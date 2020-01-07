#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o OR4D6.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i OR4D6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR4D9 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR4E2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR4F13P 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR4F15 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR4F17 
