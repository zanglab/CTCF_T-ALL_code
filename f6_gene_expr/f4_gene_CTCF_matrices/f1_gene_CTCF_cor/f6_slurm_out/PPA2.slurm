#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PPA2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PPA2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PPAN 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PPAN-P2RY11 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PPARA 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PPARD 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PPARG 
