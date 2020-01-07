#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o EGLN2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i EGLN2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EGLN3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EGR1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EGR2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EGR3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EGR4 
