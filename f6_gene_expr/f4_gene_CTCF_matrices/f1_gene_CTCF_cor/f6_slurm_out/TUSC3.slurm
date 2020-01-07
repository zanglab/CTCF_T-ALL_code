#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TUSC3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TUSC3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TUSC5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TUT1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TVP23A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TVP23B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TVP23C 
