#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SLAMF1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SLAMF1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLAMF6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLAMF7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLAMF8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLAMF9 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLBP 
