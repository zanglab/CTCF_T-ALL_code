#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o IRF9.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i IRF9 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IRGC 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IRGM 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IRGQ 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IRS1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IRS2 
