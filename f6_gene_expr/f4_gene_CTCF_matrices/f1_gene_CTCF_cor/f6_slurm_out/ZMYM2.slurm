#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ZMYM2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ZMYM2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZMYM3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZMYM4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZMYM5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZMYM6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZMYND10 
