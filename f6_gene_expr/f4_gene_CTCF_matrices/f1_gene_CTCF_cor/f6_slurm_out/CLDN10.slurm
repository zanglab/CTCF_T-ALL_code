#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CLDN10.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CLDN10 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CLDN11 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CLDN12 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CLDN14 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CLDN15 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CLDN16 
