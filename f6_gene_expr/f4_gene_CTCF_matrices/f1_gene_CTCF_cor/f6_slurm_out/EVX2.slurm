#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o EVX2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i EVX2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EWSR1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EXD1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EXD2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EXD3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EXO1 
