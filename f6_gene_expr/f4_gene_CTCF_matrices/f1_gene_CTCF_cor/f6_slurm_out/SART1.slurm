#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SART1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SART1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SART3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SASH1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SASH3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SASS6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SAT1 
