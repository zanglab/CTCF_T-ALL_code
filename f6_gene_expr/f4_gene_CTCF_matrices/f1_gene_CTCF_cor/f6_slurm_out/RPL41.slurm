#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o RPL41.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i RPL41 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RPL5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RPL6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RPL7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RPL7A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RPL7L1 
