#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o FOXE1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i FOXE1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FOXE3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FOXF1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FOXF2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FOXG1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FOXH1 
