#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o NFASC.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i NFASC 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NFAT5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NFATC1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NFATC2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NFATC2IP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NFATC3 
