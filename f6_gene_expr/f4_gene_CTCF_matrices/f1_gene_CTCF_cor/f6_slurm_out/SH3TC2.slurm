#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SH3TC2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SH3TC2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SH3YL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SHANK1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SHANK2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SHANK3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SHARPIN 
