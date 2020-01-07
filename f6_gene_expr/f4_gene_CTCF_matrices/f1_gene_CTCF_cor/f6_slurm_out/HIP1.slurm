#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o HIP1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i HIP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HIP1R 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HIPK1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HIPK2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HIPK3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HIPK4 
