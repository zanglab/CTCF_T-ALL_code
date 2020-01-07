#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o MPP4.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i MPP4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MPP5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MPP6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MPP7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MPPE1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MPPED1 
