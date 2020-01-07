#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TIAM1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TIAM1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TIAM2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TICAM1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TICAM2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TICRR 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TIE1 
