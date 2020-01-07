#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ZDHHC22.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ZDHHC22 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZDHHC23 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZDHHC24 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZDHHC3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZDHHC4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZDHHC5 
