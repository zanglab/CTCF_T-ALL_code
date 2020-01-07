#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CLDN23.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CLDN23 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CLDN24 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CLDN3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CLDN34 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CLDN4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CLDN5 
