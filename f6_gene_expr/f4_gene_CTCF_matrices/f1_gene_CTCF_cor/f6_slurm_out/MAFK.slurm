#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o MAFK.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i MAFK 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MAG 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MAGEA1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MAGEA10 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MAGEA11 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MAGEA12 
