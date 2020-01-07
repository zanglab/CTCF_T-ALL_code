#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TUBGCP4.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TUBGCP4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TUBGCP5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TUBGCP6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TUFM 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TUFT1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TULP1 
