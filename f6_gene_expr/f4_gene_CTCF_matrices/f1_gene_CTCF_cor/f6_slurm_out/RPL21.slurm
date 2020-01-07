#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o RPL21.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i RPL21 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RPL21P28 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RPL21P44 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RPL22 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RPL22L1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RPL23 
