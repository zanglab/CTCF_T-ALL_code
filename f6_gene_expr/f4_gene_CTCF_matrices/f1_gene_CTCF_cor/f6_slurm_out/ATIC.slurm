#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ATIC.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ATIC 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ATL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ATL2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ATL3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ATM 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ATMIN 
