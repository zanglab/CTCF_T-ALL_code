#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CLDN17.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CLDN17 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CLDN18 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CLDN19 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CLDN2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CLDN20 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CLDN22 
