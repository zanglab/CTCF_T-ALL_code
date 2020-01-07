#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o OR4K2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i OR4K2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR4M1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR4M2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR4N2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR4N3P 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR4N5 
