#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ACSBG1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ACSBG1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ACSBG2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ACSF2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ACSF3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ACSL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ACSL3 
