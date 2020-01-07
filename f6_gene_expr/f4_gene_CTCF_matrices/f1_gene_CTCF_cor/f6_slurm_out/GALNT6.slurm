#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o GALNT6.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i GALNT6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GALNT7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GALNT8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GALNT9 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GALNTL5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GALNTL6 
