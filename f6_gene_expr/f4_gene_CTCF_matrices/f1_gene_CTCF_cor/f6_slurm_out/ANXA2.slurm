#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ANXA2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ANXA2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ANXA2P1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ANXA2P2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ANXA2P3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ANXA2R 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ANXA3 
