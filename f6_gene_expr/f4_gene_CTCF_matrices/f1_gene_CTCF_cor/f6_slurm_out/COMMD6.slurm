#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o COMMD6.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i COMMD6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i COMMD7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i COMMD8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i COMMD9 
time python ../7_cor_CTCF_binding_GeneExpr.py -i COMP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i COMT 
