#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o FRG1HP.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i FRG1HP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FRG1JP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FRG2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FRG2B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FRG2C 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FRG2DP 
