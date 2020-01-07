#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o MAGOHB.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i MAGOHB 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MAGT1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MAK 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MAK16 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MAL 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MAL2 
