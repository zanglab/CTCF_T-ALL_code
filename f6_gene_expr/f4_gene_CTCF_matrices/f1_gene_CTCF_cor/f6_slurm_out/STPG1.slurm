#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o STPG1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i STPG1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i STPG2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i STRA6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i STRA8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i STRADA 
time python ../7_cor_CTCF_binding_GeneExpr.py -i STRADB 
