#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ABCB1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ABCB1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ABCB10 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ABCB11 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ABCB4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ABCB5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ABCB6 
