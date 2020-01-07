#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o IQCD.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i IQCD 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IQCE 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IQCF1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IQCF2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IQCF3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IQCF5 
