#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CCP110.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CCP110 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CCPG1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CCR1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CCR10 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CCR2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CCR3 
