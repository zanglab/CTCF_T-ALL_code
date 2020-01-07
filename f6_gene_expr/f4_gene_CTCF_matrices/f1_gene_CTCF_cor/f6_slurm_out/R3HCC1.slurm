#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o R3HCC1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i R3HCC1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i R3HCC1L 
time python ../7_cor_CTCF_binding_GeneExpr.py -i R3HDM1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i R3HDM2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i R3HDM4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i R3HDML 
