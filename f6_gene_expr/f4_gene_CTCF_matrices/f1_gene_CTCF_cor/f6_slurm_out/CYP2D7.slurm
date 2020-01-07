#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CYP2D7.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CYP2D7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CYP2E1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CYP2F1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CYP2G1P 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CYP2J2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CYP2R1 
