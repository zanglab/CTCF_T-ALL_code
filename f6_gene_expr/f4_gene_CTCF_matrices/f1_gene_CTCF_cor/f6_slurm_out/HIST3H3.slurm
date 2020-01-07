#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o HIST3H3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i HIST3H3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HIST4H4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HIVEP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HIVEP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HIVEP3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HJURP 
