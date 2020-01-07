#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ISG20L2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ISG20L2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ISL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ISL2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ISLR 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ISLR2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ISM1 
