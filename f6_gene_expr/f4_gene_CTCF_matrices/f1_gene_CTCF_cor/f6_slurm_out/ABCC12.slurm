#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ABCC12.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ABCC12 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ABCC13 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ABCC2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ABCC3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ABCC4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ABCC5 
