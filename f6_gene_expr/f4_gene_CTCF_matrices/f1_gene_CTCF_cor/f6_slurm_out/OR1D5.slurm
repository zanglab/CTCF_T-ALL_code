#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o OR1D5.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i OR1D5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR1E1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR1E2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR1F1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR1F2P 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR1G1 
