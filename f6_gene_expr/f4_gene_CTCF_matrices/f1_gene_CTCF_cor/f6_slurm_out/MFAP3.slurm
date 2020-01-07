#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o MFAP3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i MFAP3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MFAP3L 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MFAP4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MFAP5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MFF 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MFGE8 
