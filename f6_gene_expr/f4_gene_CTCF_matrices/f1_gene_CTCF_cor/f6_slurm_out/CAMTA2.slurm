#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CAMTA2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CAMTA2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CAND1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CAND2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CANT1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CANX 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CAP1 
