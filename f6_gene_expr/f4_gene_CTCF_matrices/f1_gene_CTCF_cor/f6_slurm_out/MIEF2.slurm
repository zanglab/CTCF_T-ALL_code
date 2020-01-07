#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o MIEF2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i MIEF2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MIEN1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MIER1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MIER2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MIER3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MIF 
