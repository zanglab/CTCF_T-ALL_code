#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TYRO3P.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TYRO3P 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TYROBP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TYRP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TYSND1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TYW1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TYW1B 
