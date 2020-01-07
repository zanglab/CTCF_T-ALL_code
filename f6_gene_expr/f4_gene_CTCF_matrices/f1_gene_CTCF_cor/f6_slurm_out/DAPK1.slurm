#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o DAPK1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i DAPK1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DAPK2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DAPK3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DAPL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DAPP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DARS 
