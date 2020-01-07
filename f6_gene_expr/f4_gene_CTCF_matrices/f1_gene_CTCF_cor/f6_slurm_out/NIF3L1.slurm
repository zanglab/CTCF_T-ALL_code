#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o NIF3L1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i NIF3L1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NIFK 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NIM1K 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NIN 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NINJ1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NINJ2 
