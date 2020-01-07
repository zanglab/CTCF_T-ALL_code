#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PHC1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PHC1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PHC2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PHC3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PHEX 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PHF1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PHF10 
