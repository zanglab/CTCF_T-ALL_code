#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o RPL23AP87.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i RPL23AP87 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RPL23P8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RPL24 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RPL26 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RPL26L1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RPL27 
