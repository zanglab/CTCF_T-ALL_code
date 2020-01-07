#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CAMK2N2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CAMK2N2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CAMK4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CAMKK1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CAMKK2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CAMKMT 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CAMKV 
