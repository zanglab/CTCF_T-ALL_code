#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o THAP2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i THAP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i THAP3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i THAP4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i THAP5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i THAP6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i THAP7 
