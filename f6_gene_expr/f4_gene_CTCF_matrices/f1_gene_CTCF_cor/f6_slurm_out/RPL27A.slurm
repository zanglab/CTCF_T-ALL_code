#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o RPL27A.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i RPL27A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RPL28 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RPL29 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RPL29P2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RPL3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RPL30 
