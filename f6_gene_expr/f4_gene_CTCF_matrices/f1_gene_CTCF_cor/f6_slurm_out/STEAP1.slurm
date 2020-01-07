#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o STEAP1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i STEAP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i STEAP1B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i STEAP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i STEAP3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i STEAP4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i STH 
