#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SCGB3A2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SCGB3A2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SCGN 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SCHIP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SCIMP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SCIN 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SCLT1 
