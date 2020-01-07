#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ME2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ME2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ME3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MEA1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MEAF6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MECOM 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MECP2 
