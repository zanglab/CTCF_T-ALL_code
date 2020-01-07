#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o BGN.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i BGN 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BHLHA15 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BHLHA9 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BHLHB9 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BHLHE22 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BHLHE23 
