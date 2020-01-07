#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o RGS8.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i RGS8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RGS9 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RGS9BP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RGSL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RHAG 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RHBDD1 
