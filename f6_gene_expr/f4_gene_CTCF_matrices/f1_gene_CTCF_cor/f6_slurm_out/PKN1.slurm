#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PKN1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PKN1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PKN2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PKN3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PKNOX1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PKNOX2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PKP1 
