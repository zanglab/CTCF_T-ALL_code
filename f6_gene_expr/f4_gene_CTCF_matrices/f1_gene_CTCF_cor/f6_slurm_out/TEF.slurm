#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TEF.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TEF 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TEFM 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TEK 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TEKT1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TEKT2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TEKT3 
