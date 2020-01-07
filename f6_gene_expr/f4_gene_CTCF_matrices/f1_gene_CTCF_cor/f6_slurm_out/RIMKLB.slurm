#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o RIMKLB.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i RIMKLB 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RIMS1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RIMS2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RIMS3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RIMS4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RIN1 
