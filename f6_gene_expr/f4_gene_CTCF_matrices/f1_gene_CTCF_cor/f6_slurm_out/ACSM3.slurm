#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ACSM3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ACSM3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ACSM4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ACSM5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ACSM6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ACSS1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ACSS2 
