#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TCP11L1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TCP11L1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TCP11L2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TCP11X2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TCTA 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TCTE1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TCTE3 
