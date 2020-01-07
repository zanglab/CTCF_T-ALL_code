#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ADAMTS5.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ADAMTS5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ADAMTS6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ADAMTS7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ADAMTS7P1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ADAMTS8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ADAMTS9 
