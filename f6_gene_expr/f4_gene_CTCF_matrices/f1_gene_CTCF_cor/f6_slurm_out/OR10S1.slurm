#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o OR10S1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i OR10S1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR10T2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR10V1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR10V2P 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR10W1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR10X1 
