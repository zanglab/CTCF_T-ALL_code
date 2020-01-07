#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o XKRY.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i XKRY 
time python ../7_cor_CTCF_binding_GeneExpr.py -i XKRY2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i XPA 
time python ../7_cor_CTCF_binding_GeneExpr.py -i XPC 
time python ../7_cor_CTCF_binding_GeneExpr.py -i XPNPEP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i XPNPEP2 
