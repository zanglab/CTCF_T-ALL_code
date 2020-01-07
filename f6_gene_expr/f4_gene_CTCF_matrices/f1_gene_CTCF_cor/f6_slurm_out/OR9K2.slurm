#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o OR9K2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i OR9K2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR9Q1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR9Q2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ORAI1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ORAI2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ORAI3 
