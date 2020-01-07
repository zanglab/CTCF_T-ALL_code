#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TDG.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TDG 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TDGF1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TDGF1P3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TDH 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TDO2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TDP1 
