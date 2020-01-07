#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o OR1L3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i OR1L3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR1L4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR1L6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR1L8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR1M1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR1N1 
