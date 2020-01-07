#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SHISA2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SHISA2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SHISA3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SHISA4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SHISA5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SHISA6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SHISA7 
