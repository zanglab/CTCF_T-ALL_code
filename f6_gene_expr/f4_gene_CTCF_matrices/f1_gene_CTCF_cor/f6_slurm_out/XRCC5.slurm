#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o XRCC5.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i XRCC5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i XRCC6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i XRCC6P5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i XRN1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i XRN2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i XRRA1 
