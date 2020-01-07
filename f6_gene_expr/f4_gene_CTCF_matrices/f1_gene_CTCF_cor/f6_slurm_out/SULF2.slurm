#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SULF2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SULF2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SULT1A1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SULT1A2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SULT1A3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SULT1A4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SULT1B1 
