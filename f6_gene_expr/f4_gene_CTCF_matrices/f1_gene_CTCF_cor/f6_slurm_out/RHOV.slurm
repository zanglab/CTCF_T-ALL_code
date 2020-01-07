#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o RHOV.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i RHOV 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RHOXF1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RHOXF1P1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RHOXF2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RHOXF2B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RHPN1 
