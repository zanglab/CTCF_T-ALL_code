#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o RECQL5.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i RECQL5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i REEP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i REEP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i REEP3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i REEP4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i REEP5 
