#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o OR51F2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i OR51F2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR51G1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR51G2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR51I1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR51I2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR51L1 
